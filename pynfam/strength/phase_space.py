# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import zip
from builtins   import range
from builtins   import object
# -------------- Utilities -----------------
from scipy.special import loggamma, factorial2
import copy
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import barycentric_interpolate, BarycentricInterpolator
# ----------- Relative Imports -------------
from ..config import ALPHA, HBARC, HBAR_MEC, MEC2, KAPPA
from ..config import DMASS_NH, R0, MN, EPSILON
from ..config import DEFAULTS
from .electron_capture import kappa_x, betasq_x, energy_x, chem_pot
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#=========================================================================#
#                           Class phaseSpace                              #
#=========================================================================#
class phaseSpace(object):
    """
    An instance of the class phaseSpace contains the methods for computing
    phase space integrals and fits to those integrals which can be analytically
    continued to the complex plane.

    This is designed to work similarly to the fortProcess objects, in that
    it has default settings that are overridden.

    Args:
        beta (str): The beta decay type (ONLY '-' IS IMPLEMENTED)

    Attributes:
        beta (str): The beta decay type.
        _settings (dict): The underlying phase space settings.
    """

    def __init__(self, beta):
        # type of beta decay ('+', '-', 'c')
        self.beta = beta
        if beta not in [u'-',u'+',u'c']:
            raise ValueError(u"Invalid beta decay type requested in phaseSpace.")

        # default settings
        self._settings = self._getDefaults()

    @property
    def approx(self):
        """ approx (str): The type of approximation used to analytically continue
        the phase space integrals.
        """
        return self._settings[u'psi_approx']

    @property
    def GA(self):
        """ GA (float): Axial vector coupling, gA ~ -1.27
        """
        return self._settings[u'GA']

    @property
    def GV(self):
        """ GV (float): Vector coupling. gV = +1.0 but can be overridden
        """
        return self._settings[u'GV']

    @property
    def LAM(self):
        """ LAM (float): Positive parameter -gA/gV entering shape factor
        from 'Electron Radial Wavefunctions', 1982, Behrens and Buhring
        """
        return abs(self._settings[u'GA']/self._settings[u'GV'])


    #-----------------------------------------------------------------------
    def updateSettings(self, override):
        """
        Update the settings for the phase space calculations.

        Since psi_approx is itself available as override setting, and each value
        of psi_approx has different settings (like the famContour) we have to
        update the available settings if we change it.

        Args:
            override (dict): Settings to override.
        """
        #if u'psi_approx' in list(override.keys()):
        #    self._settings = self._getDefaults(override[u'psi_approx'])

        for h in override:
            if h not in list(self._settings.keys()):
                raise KeyError(u"Invalid override setting '{:}' ".format(h) +\
                        u"for phase space approx '{:}'.".format(self.approx))
            else:
                self._settings[h] = override[h]


    #-----------------------------------------------------------------------
    def analyticPsi(self, fpsi, w0min, w0max, debug=False, mu=None):
        """
        Constructs an analytic approximation to phase space integrals on the real axis
        given a function which computes the integrals numerically.

        Two types of approximations are implemented:

            * RATINT: A rational function interpolation of the integrals on the real axis.
            * POLYFIT: A polynomial least squares fit of the integrals on the real axis.

        Both types approximate the phase space integrals as a function of the maximum
        beta particle energy W0, which is itself a simple function of the QRPA energy.
        The fitting data are computed on a chebychev grid to reduce fit instabilities
        (e.g. Runge phenomenon), and the phase space integrals are computed with
        Gauss-Legendre quadrature.
        """

        # Get chebychev grid of 'x-values' (w0) along real axis for the approximation
        if self.approx == u'RATINT':
            npts_approx = self._settings[u'ratint_pts']
        elif self.approx == u'POLYFIT':
            npts_approx = self._settings[u'polynomial_fit_pts']

        w0_vals, cwts = np.polynomial.chebyshev.chebgauss(npts_approx)
        w0_vals, cwts = transform_interval(w0_vals,cwts,w0min,w0max)

        # Calculate the corresponding 'y-values' (phase space integrals)
        psi = fpsi(w0_vals)

        # Contruct the approximation using known x=w0grid and y=f_n(w0)
        if self.approx == u'RATINT':
            psi_approx = thieleInterpolator(w0_vals, psi)
        elif self.approx == u'POLYFIT':
            order = self._settings[u'polynomial_order']
            if mu is None:
                psi_approx = polynomialFit(w0_vals, psi, order)
            else:
                psi_approx = polyexpFit(w0_vals, psi, order, mu)

        # Return the wrapper fct (+ grid used for approx if debug on)
        if not debug:
            return psi_approx
        else:
            return psi_approx, w0_vals, psi


    #-----------------------------------------------------------------------
    def psiFct(self, n, Zd, A, eqrpamax, eqrpamin, approx=True, Evloss=False):
        """
        Construct an analytic approximation of the phase space for BETA DECAY.
        Charge screening is applied for Zd<0, which should be supplied in the case of beta-plus.
        My tests show that screening is not important for beta-minus.

        Args:
            n (int): The index for the phase space type (Mustonent et al., Phys. Rev. C 90, 024308 (2014)).
            Zd (int): Charge of the daughter nucleus.
            A (int): Mass of the daughter/parent nucleus.
            eqrpamax (float): Maximum QRPA energy for beta decay.
            eqrpamin (float): Minimum QRPA energy for phase space integral approximation (default 0)
            debug (bool): If true, also returns the x=w0, y=psi values used for the approximation.

        Returns:
            ufunc
            ndarray: The x=w0 values used for the approximation (only if debug=True).
            ndarray: The y=psi values used for the approximation (only if debug=True).
        """
        # Force the correct inputs for beta plus and ec
        if self.beta == u'-':
            Zd = abs(Zd)
            sc = False
        elif self.beta == u'+':
            Zd = -abs(Zd)
            sc = True
        else:
            raise RuntimeError("psiFct called for beta type={:} not in [-,+].".format(self.beta))

        # Override screening if requested
        sc_user = self._settings[u'screening']
        if sc_user is not None: sc = sc_user

        # Max beta particle energies
        w0min = 1.0 # At EQRPAmax, just enough to create the particles mass
        w0max = 1.0 + ((eqrpamax-eqrpamin)/MEC2) # Interpolate over whole interval, possibly to EQRPA<0.

        # Construct the approximation
        def fpsi(w0): return self.calcPsi(n, Zd, A, w0, sc=sc, debug=False, Evloss=Evloss)
        if approx:
            psi_fct = self.analyticPsi(fpsi, w0min, w0max)
        else:
            psi_fct = fpsi # Note, real argument (w0) is enforced

        return psi_fct


    #-----------------------------------------------------------------------
    def psecFct(self, n, Zp, A, Bx_shift=False, eqrpamax=None, interval=None, temper=None, approx=True, Evloss=False):
        """
        Construct an analytic approximation of the phase space for ELECTRON CAPTURE.

        Terrestrial EC is already a polynomial in w0 (thus no approximation needed),
        but finite temp EC is an integral function like for beta decay.
        """

        if self.beta == u'c':
            Zp = abs(Zp) # parent
            Zd = Zp - 1  # daughter
            Bk = 0.0     # binding energy shift
        else:
            raise RuntimeError("psiFct called for beta type={:} != 'c'.".format(self.beta))

        if temper is None:
            # EQRPAmax neglects electron binding energy Bk.
            # w0 = (eqrpamax_ec - eqrpa)/mec2 - 1   and    pv = w0 + 1 - Bk
            # Optionally shift w0 to use the correct EQRPAmax' = EQRPAmax - Bk
            # For contour integrals, the rate incurs a small error on the order of
            # Bk^2 or Bk^4 times the strength between EQRPAmax' and EQRPAmax.
            # Bk is approximately <~ 0.2 MeV for most nuclei, so this error is small.
            if Bx_shift:
                Zk = Zd - 0.35 # Effective daughter charge due to screening
                Bk = 1 - np.sqrt(1 - (ALPHA*Zk)**2) # Hydrogen-like K-shell BE

            # EC phase space is a polynomial of w0 so is already analytic
            def psec_fct(w0in): return self.calcPsec(n, Zp, A, w0in-Bk)
        else:
            # Override screening if requested
            sc = False
            sc_user = self._settings[u'screening']
            if sc_user is not None: sc = sc_user

            # Phase space integrals as a function
            logpYe = self._settings[u'log(pYe)']
            logpYe = np.array(logpYe, copy=False, ndmin=1)
            if logpYe.size > 1:
                raise RuntimeError("psecFct got multiple densities for FTEC.")
            mu = chem_pot(logpYe, temper)[0] # mu[MeV]
            def fpsi(w0): return self.calcPsecft(n, Zp, A, w0, temper, mu, sc=sc, Evloss=Evloss)

            # Construct the approximation
            if approx:
                w0min = (eqrpamax - interval[1])/MEC2 - 1.0
                w0max = (eqrpamax - interval[0])/MEC2 - 1.0
                psec_fct = self.analyticPsi(fpsi, w0min, w0max, mu=mu)
            else:
                psec_fct = fpsi

        return psec_fct


    #-----------------------------------------------------------------------
    @staticmethod
    def psi_integrand(wx, n, Zd, A, w0, sc=False, mu_T=None, noFermi=False, useG=False, Evloss=False):
        """ The standard beta decay phase space integrand used for beta-plus/beta-minus/FT-EC. """

        wx = np.asarray(wx)
        f = np.zeros_like(wx, dtype=np.float_)

        # Avoid invalid values (for 0T open contours above EQRPAmax)
        wx[wx<1.0] = 1.0

        px = np.sqrt(wx ** 2.0 - 1.0)
        if mu_T is None:
            # Beta+ or Beta-
            pv = w0 - wx
            fe = 1
        else:
            # FTEC
            mu, temp = mu_T
            pv = w0 + wx
            fe = 0.5*(1.0-np.tanh(0.5*((wx-1)*MEC2 - mu)/temp)) # (wx-1)*mec2 = Electron KE in MeV

        if n == 1:
            gn = gam_ke(1, Zd) * mu_ke(1)
        elif n == 2:
            gn = wx
        elif n == 3:
            gn = wx * wx
        elif n == 4:
            gn = wx * wx * wx
        elif n == 5:
            gn = wx * pv * pv
        elif n == 6:
            if noFermi:
                # lamk = c*F{k-1}/F0,  lamk_F0Fk == lamk*F0/F{k-1} = c
                # ps ~ lam2 F0 L0 = (lam2 F0/F1) F1 L0 = c F1 L0, then noFermi = remove F1 L0
                lke = lambda_ke_F0Fk(2, Zd)
            else:
                lke = lambda_ke(2, Zd, A, wx, sc)
            gn = wx * px * px * lke

        if noFermi:
            F = 1
        else:
            F = Fermi(0, Zd, A, wx, sc) *  L0(Zd)
            if useG:
                F = F * px / wx # This is now G

        # Neutrino energy loss is same calculation but with pv^3
        if Evloss:
            fpv = pv*pv*pv
        else:
            fpv = pv*pv

        if not useG:
            f = px * fpv * gn * fe * F # the usual
        else:
            f = wx * fpv * gn * fe * F # sub wG = pF (G=p/w F varies slowly, better behaved)

        f[np.logical_not(np.isfinite(f))] = 0

        return f


    #-----------------------------------------------------------------------
    def calcPsi(self, n, Zd, A, w0, sc=False, debug=False, Evloss=False):
        """
        Calculate the phase space integral from 1 to w0 with Gauss-
        Legendre quadrature on psi_glpts sized grid (from settings).

        Args:
            n (int): The index for the phase space type (Mustonent et al., Phys. Rev. C 90, 024308 (2014)).
            Zd (int): Charge of the daughter nucleus.
            A (int): Mass of the daughter/parent nucleus.
            w0 (float): Maximum electron energy for the decay :math:`W_0=1+(EQRPA_{max} - EQRPA)/m_ec^2`
            debug (bool): If True, also return the x=w and y=phase space values used for the integration.

        Returns:
            float
            ndarray: The x=w values used for the integration (only if debug=True).
            ndarray: The y=phase space used for the integration (only if debug=True).
        """

        w0 = np.asarray(w0)

        if n not in list(range(1,7)):
            raise ValueError(u"Phase space factor n takes values 1 to 6.")
        # PSI is not complex analytic so should receive real valued arguments
        if not np.isreal(w0).all():
            raise ValueError("calcPsi got complex argument w0.")

        glpts = self._settings[u'psi_glpts']
        wmin = 1.0
        def integrate(w0):
            w, glwts = np.polynomial.legendre.leggauss(glpts)
            w, glwts = transform_interval(w, glwts, wmin, w0)
            farr = self.psi_integrand(w, n, Zd, A, w0, sc=sc, Evloss=Evloss)
            fint = np.sum(farr*glwts)
            if not debug:
                return fint
            else:
                return fint, w, farr

        # Allow passing arrays of w0
        return np.vectorize(integrate)(w0)


    #-----------------------------------------------------------------------
    def calcPsec(self, n, Zp, A, w0):
        """ Calculate terrestrial EC phase space. It is just a constant depending
            on the shell, rather than an integral.
        """

        w0 = np.asarray(w0)
        gn = np.ones_like(w0, dtype=np.float_)
        # Allow complex values, but force int to float
        if np.iscomplex(w0).any():
            f = np.zeros_like(w0, dtype=np.complex_)
        else:
            f = np.zeros_like(w0, dtype=np.float_)

        if n not in list(range(1,7)):
            raise ValueError(u"Phase space factor n takes values 1 to 6.")

        # K-shell quantities
        Bd = energy_x(1,0,1,Zp-1)
        Bp = energy_x(1,0,1,Zp)
        kap = kappa_x(1,0,1)
        betax = betasq_x(1,0,1,Zp-1,A)
        nx = 1.0

        # Electron and neutrino quantities
        wx = 1.0 - Bp # Fixed electron energy
        # Behrens says W^2-1, but this would be imaginary....
        # Bambynek says 1-W^2, but this doesn't make sense...
        px = np.sqrt(1.0 - wx**2.0) # Fixed electron momentum
        pv = w0 + 1 # Note w0 = w0-Ex if Bx_shift, so pv = w0 + 1 - Ex

        # To align with beta-plus/beta-minus notation
        # (careful to use correct multiplication casting rule to accept type float or complex wx)
        if n == 1:
            gn = (gn * -np.sign(kap))
        elif n==2:
            gn = (gn * 1)
        elif n==3:
            gn = (gn * wx)
        elif n==4:
            gn = (gn * (wx * wx))
        elif n==5:
            gn = (gn * (pv * pv))
        elif n==6:
            gn = (gn * (px * px))

        # psEC is complex analytic but should technically be forced to zero at
        # unphysical neutrino momenta (pv<0). Enforce this on the real axis but
        # not for complex pv (this then incurs an error if real(pv)<0)
        if np.isreal(w0).all():
            mask = np.logical_not(w0 < -1.0)
        else:
            mask = np.ones_like(w0, dtype=np.bool) # use whole array

        f[mask] = gn[mask] * pv[mask] * pv[mask] * betax * (np.pi*0.5) * nx

        return f


    #-----------------------------------------------------------------------
    def calcPsecft(self, n, Zd, A, w0, temp, mu, sc=False, debug=False, Evloss=False):
        """
        Calculate the FT EC phase space integral from wth to infinity with
        Gauss-Legendre + Gauss-Laguerre quadrature on psi_glpts sized grid (each).

        Args:
            n (int): The index for the phase space type (Mustonent et al., Phys. Rev. C 90, 024308 (2014)).
            Zd (int): Charge of the parent nucleus.
            A (int): Mass of the daughter/parent nucleus.
            w0 (float): Maximum positron energy for the decay :math:`W_0=(EQRPA_{max-ec} - EQRPA)/m_ec^2 - 1`
            temp (float): Temperature in MeV
            mu (float): Stellar electron chemical potential (without electron rest mass, and unitless (mu=Uf/mec^2)
            sc (bool): Use charge screening
            debug (bool): If True, also return the x=w and y=phase space values used for the integration.

        Returns:
            float
            ndarray: The x=w values used for the integration (only if debug=True).
            ndarray: The y=phase space used for the integration (only if debug=True).
        """

        w0 = np.asarray(w0)
        mu = float(mu)
        tol = 0.1 # This helps the gauss lageurre capture most of the peak for small temps when the drop off is steep

        if n not in list(range(1, 7)):
            raise ValueError(u"Phase space factor n takes values 1 to 6.")
        if not np.isreal(w0).all():
            raise ValueError("calcPsecft got complex argument w0.")

        glpts = self._settings[u'psi_glpts']
        def integrate(w0):
            # Threshold is higher if neutrino momentum < 0.
            # Set the lower bound as threshold energy, upper bound as fermi energy
            if w0 < -1:
                S0 = abs(w0)
            else:
                S0 = 1
            S1 = mu/MEC2 + 1 + tol # Note: Min(mu/MEC2) = -1

            if S1 <= S0:
                # If the threshold is greater than the fermi energy, only Gauss-laguerre
                farr1 = np.zeros(1)
                w1    = np.zeros(1)
                fint1 = 0.0
                S1 = S0
            else:
                # Integrate the 1st portion (up to chem pot) with Gauss-legendre
                w1, glwts1 = np.polynomial.legendre.leggauss(glpts)
                w1, glwts1 = transform_interval(w1, glwts1, S0, S1)
                farr1 = self.psi_integrand(w1, n, Zd, A, w0, sc=sc, mu_T=(mu, temp), noFermi=False, useG=True, Evloss=Evloss)
                fint1 = np.sum(farr1*glwts1, axis=0)

            # Integrate the 2nd portion with Gauss-laguerre
            w2, glwts2 = np.polynomial.laguerre.laggauss(glpts)
            glwts2 = glwts2*np.exp(w2)
            w2 = w2 + S1
            farr2 = self.psi_integrand(w2, n, Zd, A, w0, sc=sc, mu_T=(mu, temp), noFermi=False, useG=True, Evloss=Evloss)
            fint2 = np.sum(farr2*glwts2, axis=0)

            fint = fint1 + fint2

            if debug:
                return fint, w1, farr1, w2, farr2
            else:
                return fint

        # Allow passing arrays of w0
        return np.vectorize(integrate)(w0)


    #-----------------------------------------------------------------------
    def get_max_eqrpa_ecft(self, temper, Zp, A, eqrpamax, min_eqrpa, max_eqrpa, step=0.1, tol=1e-10):
        """ Find w0 where phase space becomes zero due to machine precision. Root finding
        is not stable because the function rapidly hits zero at large negative w0s
        so brute force it. """

        # Override screening if requested
        sc = False
        sc_user = self._settings[u'screening']
        if sc_user is not None: sc = sc_user

        # Chemical potential
        logpYe = self._settings[u'log(pYe)']
        mus = chem_pot(logpYe, temper) # mu[MeV], array of len(logpYe)

        # Find the lagest EQRPA that gives a non-zero phase space
        # (among all 6 phase space types, this will be the slowest varying one, f1)
        emax = []
        for i, mu in enumerate(mus):
            eqrpa = max_eqrpa
            while eqrpa > min_eqrpa:
                w0 = (eqrpamax - eqrpa)/MEC2 - 1.0
                psi = self.calcPsecft(1, Zp, A, w0, temper, mu, sc=sc)
                if psi > tol:
                    break
                eqrpa -= step
            # Handle endpoints correctly
            if eqrpa < min_eqrpa:
                eqrpa = min_eqrpa

            # Don't get too close to FT pole
            if abs(eqrpa) < step:
                eqrpa = 0.0
            emax.append(eqrpa)

        return np.array(emax), mus


    #-----------------------------------------------------------------------
    def _getDefaults(self):#, psi_approx=None):
        """
        Return the dictionary for default phase space settings.

        Args:
            psi_approx (str): The default approximation type (default None --> config.py value).

        Returns:
            dict
        """

        default_raw = copy.deepcopy(DEFAULTS[u'psi'])

        default_settings = {}
        for h in default_raw:
            default_settings.update(default_raw[h])

        #default_settings = default_raw[u'PSI']
        #if psi_approx is None:
        #    psi_approx = default_settings[u'psi_approx']
        #default_settings.update(default_raw[psi_approx])
        #default_settings[u'psi_approx'] = psi_approx # this is redundant, but just to be safe

        return default_settings


#=========================================================================#
#                      Numerical functions                                #
#=========================================================================#
def transform_interval(x, wts, xmin, xmax):
    """
    Given gauss legendre grid and weights on -1 to 1, transform
    to interval xmin to xmax.

    Args:
        x (ndarray): The Gauss-Legendre grid on -1 to 1.
        wts (ndarray): The Gauss-Legendre weights on -1 to 1 grid.
        xmin (float): New lower bound.
        xmax (float): New upper bound.

    Returns:
        ndarray: Transformed Gauss-Legendre grid.
        ndarray: Transformed Gauss-Legendre weights.
    """
    f1  = 0.5*(xmin+xmax)
    f2  = 0.5*(xmax-xmin)
    wts = wts*f2
    x   = x*f2 + f1
    return x, wts

#-----------------------------------------------------------------------
def polynomialFit(x, y, d):
    """
    Returns a function which evaluates a polynomial least squares fit
    built from x and y.

    Args:
        x (ndarray): The x values for the fit.
        y (ndarray): The y values for the fit.
        d (int): The degree of the polynomial used for the fit.

    Returns:
        ufunc
    """
    def p(xin):
        return np.polyval(np.polyfit(x,y,d), xin)
    return p

#-----------------------------------------------------------------------
def polyexpFit(x, y, order, mu):
    """
    Returns a function which evaluates a least squares fit
    built from x and y. The fit is a polynomial if the domain
    is entirely above mu, an exponential if the domain is entirely
    below mu, or a polynomial with exponential envelope if mu is
    inside the domain.

    Args:
        x (ndarray): The x values for the fit.
        y (ndarray): The y values for the fit.
        order (int): The degree of the polynomial used for the fit.
        mu (float): The fermi energy in MeV (kinetic only)

    Returns:
        ufunc
    """

    tol = 1e-16
    xtol= 0.1
    cut = -(mu/MEC2 + 1.) # w0 for the fermi energy

    lower = cut - min(x)
    upper = max(x) - cut

    # All zero
    xz = x[y>tol]
    yz = y[y>tol]
    if xz.size==0:
        return lambda xx: np.zeros_like(xx)

    # Polynomial (whole interval)
    if lower < xtol:
        fitp = np.polyfit(x, y, order)
        return lambda xx: np.polyval(fitp, xx)

    # Exponential (only y>0 can be included in fit)
    elif upper < xtol:
        fite = np.polyfit(xz, np.log(yz), 1)
        return lambda xx: np.exp(np.polyval(fite, xx))

    # Hybrid doesn't work well, fall back to poly
    else:
        fitp = np.polyfit(x, y, order)
        return lambda xx: np.polyval(fitp, xx)

    ## Hybrid (polynomial with exponential envelope)
    ## With contour splitting, this should never happen
    #else:
    #    print("WARNING: FTEC contour was not split at fermi energy, using poly-exp fit.")
    #    fitp = np.polyfit(xz, yz, order)
    #    fite = np.polyfit(xz, np.log(yz), 1)
    #
    #    if fite[1] > np.log(1e100):
    #        cexp = 1e100
    #    else:
    #        cexp = np.exp(fite[1])
    #    p0 = list(cexp*fitp)+[fite[0]]
    #
    #    def hfct_clist(xx, coefs):
    #        if len(coefs)!= len(p0): raise ValueError
    #        ps, m = coefs[:-1], coefs[-1]
    #        return np.polyval(ps, xx)*np.exp(m*xx)
    #    hfct = lambda xx, *coefs: hfct_clist(xx, coefs)
    #    try:
    #        fith, cov = curve_fit(hfct, x, y, p0=p0, maxfev=50000)
    #        return lambda xx: hfct(xx, *fith)
    #    except RuntimeError:
    #        return lambda xx: np.polyval(fitp, xx)

#-----------------------------------------------------------------------
def thieleInterpolator(x, y):
    """
    Returns a funcion which evaluates a rational function interpolation
    built from x and y using thiele continued fractions.

    Args:
        x (ndarray): The x values for the fit.
        y (ndarray): The y values for the fit.

    Returns:
        ufunc
    """
    # Calculate the continued fraction coefficients
    p = [[yi]*(len(y)-i) for i, yi in enumerate(y)]
    for i in range(len(p)-1):
        p[i][1] = (x[i] - x[i+1])/(p[i][0] - p[i+1][0] + 1e-15)
    for i in range(2, len(p)):
        for j in range(len(p)-i):
            p[j][i] = (x[j]-x[j+i])/(p[j][i-1]-p[j+1][i-1]) + p[j+1][i-2]
    p0 = p[0]

    # Evaluate the continued fraction
    def t(xin):
        a = 0
        for i in range(len(p0)-1, 1, -1):
            a = ((xin - x[i-1])/(p0[i] - p0[i-2] + a))
        return y[0] + ((xin-x[0])/(p0[1]+a))

    # Check the interpolation passes though the provided points
    for xx, yy in zip(x, y):
        diff   = abs(yy - t(xx))
        mean   = 0.5*(abs(yy + t(xx)))
        reldiff = (diff/(mean+EPSILON))
        if min(reldiff,diff) > 1e-5:
            print(u"WARNING: ratint does not pass through points precisely.")
            print(u"Relative diff: {:}, Abs diff: {:} at x={:}.".format(reldiff, diff, xx))

    # Return the evaluating function
    return t

#=========================================================================#
#                   Coulomb and Fermi Function                            #
#=========================================================================#
def NZ_tilde(Zd):
    """
    Charge screening parameter. In general it is a function of the parent
    charge, but it is slowly varying and we take it to be a constant.

    Args:
        Zd (int): Charge of daughter nucleus.

    Returns:
        float
    """
    return 1.43

#-----------------------------------------------------------------------
def V0_shift(Zd):
    """
    Positive energy shift due to charge screening from a daughter nucleus
    of charge Zd.

    Args:
        Zd (int): Charge of daughter nucleus.

    Returns:
        float
    """
    # Zd - 1 = Zp always...  B-: +|Zd| - 1 = +|Zp|,  B+: -|Zd| - 1 = -|Zp|
    # The absolute value is to avoid python taking complex root and
    # returning a complex value. Since it's really ((Zd-1)^4)^1/3 the
    # result should be invariant under sign change anyway.
    return NZ_tilde(Zd)*ALPHA**2*abs(Zd-1)**(4.0/3.0)

#-----------------------------------------------------------------------
def w_screen(Zd, w):
    r"""
    Shifted electron energy due to charge screening (denoted :math:`\tilde{w}`).

    Args:
        Zd (int): Charge of daughter nucleus.
        w (float): Electron energy.

    Returns:
        ndarray
    """
    return np.asarray(w) - V0_shift(Zd)

#-----------------------------------------------------------------------
def mu_ke(ke):
    """
    Coulomb function: Taken to be 1, only need mu_1 in shape factor.

    Args:
        ke (int): electron spherical wave expansion index.

    Returns:
        float
    """
    return 1.0

#-----------------------------------------------------------------------
def gam_ke(ke, Zd):
    """
    Coulomb function: Appears in the definition of the Fermi Functions
    and first forbiden shape factor.

    Args:
        ke (int): electron spherical wave expansion index.
        Zd (int): Charge of daughter nucleus.

    Returns:
        float
    """
    return np.sqrt(float(ke)**2.0 - (ALPHA*Zd)**2.0)

#-----------------------------------------------------------------------
def lambda_ke(ke, Zd, A, w, sc=False):
    r"""
    First approximation to the Coulomb function :math:`\lambda_k`.
    Charge screening is treated as in Suhonen NPA 563, 205 (1993).

    Args:
        ke (int): The order of the coulomb function.
        Zd (int): Charge of the daughter nucleus.
        A (int): Mass of the parent/daughter nucleus.
        w (float): The electron energy.
        sc (bool): Use charge screening

    Returns:
        float
    """

    w = np.asarray(w)

    # Initialize
    f0  = np.zeros_like(w, dtype=np.float_)
    fk  = np.zeros_like(w, dtype=np.float_)
    lke = np.zeros_like(w, dtype=np.float_)

    if not np.isreal(w).all():
        raise ValueError("lambda_ke got complex energy argument.")

    # Charge screening
    ws = w
    if sc: ws = w_screen(Zd, w)

    # Special case: when p --> 0, the shape should go to 0, so we can safely return 0
    # Also, where ws <=1, Fermi function goes to infinity, so avoid evaluation at these points
    mask = np.logical_not(np.logical_or(abs(ws-1.0) < EPSILON,  ws <= 1.0))

    f0[mask] = Fermi(0, Zd, A, w[mask], sc) # Use w here, converted to w_screen inside Fermi
    fk[mask] = Fermi(ke-1, Zd, A, w[mask], sc)
    g1 = gam_ke(1,Zd)
    gk = gam_ke(ke,Zd)
    lke[mask] = (ke+gk)/(ke*(1+g1))*(fk[mask]/f0[mask])

    return lke

#-----------------------------------------------------------------------
def lambda_ke_F0Fk(ke, Zd):
    r"""
    lambda_ke times F0/Fk, which removes the dependence on Fermi function.

    Args:
        ke (int): The order of the coulomb function.
        Zd (int): Charge of the daughter nucleus.

    Returns:
        float
    """

    g1 = gam_ke(1,Zd)
    gk = gam_ke(ke,Zd)
    return (ke+gk)/(ke*(1+g1))

#-----------------------------------------------------------------------
def L0(Zd):
    """
    Coulomb function: Simplest approximation for L0 without charge-screening.

    Args:
        Zd (int): Charge of daughter nucleus.
        sc (bool): Include corrections for charge screening.

    Returns:
        float
    """
    return 0.5*(1.0 + gam_ke(1, Zd))

#-----------------------------------------------------------------------
def L0_2(Zd, A, w, sc=False):
    r"""
    Coulomb function: Approximation to L0 for :math:`|Z| <~ 15`.

    Args:
        Zd (int): Charge of daughter nucleus.
        A (int): Mass of parent/daughter nucleus.
        w (float): Electron energy.

    Returns:
        float
    """
    w = np.asarray(w)

    if sc: raise ValueError(u"Charge screening not yet supported for L0_2.")

    R = ((R0/HBAR_MEC))*A**(1.0/3.0)
    L0_2 = 1.0 - ALPHA*Zd*w*R + (13.0/60.0)*(ALPHA*Zd)**2 - 0.5*((ALPHA*Zd*R/w))

    return L0_2

#-----------------------------------------------------------------------
def Fermi(F, Zd, A, w, sc=False):
    r"""
    Fermi Function, excluding L0 which must be multiplied separately.

    We first compute ln(F) to avoid overflow and underflow in individual
    terms. Once combined the value is reasonable to then exponentiate.

    Note:
        scipy loggamma handles complex arguments just fine, see scipy documentation
        for difference between loggamma and gammaln.

        About numerical stability:
        As :math:`w \rightarrow 1`, y become large (~10^4), thus :math:`\exp(\pi*y)` overflows.
        For large y, the Re and Im parts of :math:`\gamma(1+iy)` underflow, approaching zero.
        Taking the :math:`ln`, :math:`\ln \exp(\pi*y) \propto y` and
        :math:`\ln \gamma(1+iy) \propto -y` (not sure exact relation).
        These roughly balance when added, allowing us to exponentiate the result and
        recover the Fermi function avoiding over/underflows.

    Args:
        F (int): The type of Fermi function 0 or 1.
        Zd (int): Charge of the daughter nucleus.
        A (int): Mass of the parent/daughter nucleus.
        w (float): The electron energy.
        sc (bool): Include correction for charge screening (default False).

    Returns:
        float

    Raises:
        ValueError
    """

    w = np.asarray(w)

    # Initialize arrays
    p          = np.zeros_like(w, dtype=np.float_)
    pref       = np.zeros_like(w, dtype=np.float_)
    factor     = np.zeros_like(w, dtype=np.float_)
    p_us       = np.zeros_like(w, dtype=np.float_)
    y          = np.zeros_like(w, dtype=np.float_)
    ln_expy    = np.zeros_like(w, dtype=np.float_)
    garg       = np.zeros_like(w, dtype=np.complex_)
    ln_gamnum1 = np.zeros_like(w, dtype=np.complex_)
    ln_gamnum2 = np.zeros_like(w, dtype=np.complex_)
    ln_F       = np.zeros_like(w, dtype=np.complex_)
    result     = np.zeros_like(w, dtype=np.float_)

    if not np.isreal(w).all():
        raise ValueError("Fermi function got complex energy argument.")

    # Treat the energy w
    # ------------------------------------------------------
    # Warn if w < 1 (do not warn for screened ws=w-V0 < 1)
    # if (w < 1.0).any(): print(u"WARNING: Re(w) < 1 in Fermi.")

    # Charge screening
    if sc:
        w_us = w
        w = w_screen(Zd, w)

    # Special case: when p --> 0, shape should go to 0, so we can safely return 0
    # This also applies to screened ws=w-V0.
    mask = np.logical_not(np.logical_or(abs(w-1.0) < EPSILON,  w < 1.0))

    # Compute momentum p (after weve determined it will be real)
    # ------------------------------------------------------
    p[mask] = np.sqrt(w[mask]**2.0 - 1.0)

    # Compute Fermi function (with possible Rose screening)
    # ------------------------------------------------------
    if not sc:
        pref[mask] = 1.0
    else:
        p_us[mask] = np.sqrt(w_us[mask]**2.0 - 1.0)
        pref[mask] = ((p[mask]/p_us[mask]))*((w[mask]/w_us[mask]))

    # factor = ([ke* (2ke-1)!!]^2) * (4^ke) * (2pR)^(2gam_ke - ke)
    R = (R0/HBAR_MEC*A**(1.0/3.0))
    gF = gam_ke(F+1,Zd)
    factor[mask] = ((F+1)*factorial2((2*(F+1) - 1)))**2*4**(F+1)*(2.0*p[mask]*R)**(2.0*(gF-F-1))

    # Fermi = pre * factor * exp(pi*y) * |G(gam_ke + iy)|^2/(G(2*gam_ke + 1))^2
    y[mask] = (ALPHA*Zd*w[mask]/p[mask])
    ln_expy[mask] = np.pi*y[mask]
    garg[mask]= gF + 1j*y[mask]
    ln_gamnum1[mask] = loggamma(garg[mask])
    ln_gamnum2[mask] = np.conj(loggamma(garg[mask]))
    ln_gamden = loggamma(2.0*gF + 1.0)

    ln_F[mask]= np.log(pref[mask]) +np.log(factor[mask]) +ln_expy[mask] +ln_gamnum1[mask] +ln_gamnum2[mask] -2*ln_gamden

    # F contains abs(gamma(x+iy))^2 so should be real, thus ln F should be real
    if (np.imag(ln_F) > EPSILON).any():
        print(u"WARNING: ln(F(w)) is not exactly real, but it should be.")
    # This should never happen because 'garg' is complex and 'gF' is a sqrt and thus > 0,
    # but gamma/loggamma is undefined on the negative real axis...
    if np.isnan(ln_gamnum1).any() or np.isnan(ln_gamden).any():
        print(u"WARNING: loggamma returned NaN result in Fermi.")

    result[mask] = np.exp(np.real(ln_F[mask]))
    return result
