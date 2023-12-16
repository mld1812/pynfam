# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins import range
from builtins import object
# -------------- Utilities -----------------
import copy
import numpy as np
# ----------- Relative Imports -------------
from ..config import DEFAULTS
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS famContour                                  #
#===============================================================================#
class famContour(object):
    """
    An instance of the class famContour contains the properties and values of a
    complex energy contour on which the strength is to be computed, as well as
    the properties needed to integrate the strength along the contour.

    This is designed to work similarly to the fortProcess objects, in that
    it has default settings that are overridden.

    Args:
        contour (str): The name of the contour type.
        override (dict): Settings and values to override defaults (default None).

    Attributes:
        name (str): The name of the contour type.
        _settings (dict): Contour settings.
    """

    def __init__(self, contour, override=None):
        self.name = contour.upper()
        self.density = None
        self._settings = self._getDefaults()
        self._generateCtrData()
        if override is not None:
            self.updateSettings(override)

    @property
    def name_and_int(self):
        """ name_and_int (str): Contour name + energy interval as string
        """
        return "{:} on ({:.2e}, {:.2e})".format(self.name, self.energy_min, self.energy_max)
    @property
    def closed(self):
        """ closed (bool): Closed vs open contour.
        """
        return self._ctr_data[u'closed']
    @property
    def nr_points(self):
        """ nr_points (int): Number of discrete points along the contour.
        """
        return self._ctr_data[u'nr_points']
    @property
    def nr_compute(self):
        """ nr_compute (int): Number of points on which to compute the strength.
        """
        return self._ctr_data[u'nr_compute']
    @property
    def use_gauleg(self):
        """ use_gauleg (bool): Use a Gauss-Legendre spacing for the contour points.
        """
        return self._ctr_data[u'use_gl_ctr']
    @property
    def ctr_z(self):
        """ ctr_z (complex ndarray): The complex energy values.
        """
        return self._ctr_data[u'ctr_z']
    @property
    def ctr_dzdt(self):
        """ ctr_dzdt (complex ndarray): dz/dt where t parameterizes the contour.
        Real valued array of ones for open contours with no parameterization.
        """
        return self._ctr_data[u'ctr_dzdt']
    @property
    def theta(self):
        """ theta (ndarray): Parameterization of the contour.
        Array of zeros for open contours.
        """
        return self._ctr_data[u'theta']
    @property
    def glwts(self):
        """ glwts (ndarray): Gauss-Legendre weights.
        Array of zeros for contours not on Gauss-Legendre grid.
        """
        return self._ctr_data[u'glwts']
    @property
    def half_width(self):
        """ half_width (float): The imaginary value of the contour energies.
        None if the contour does not lie entirely at a single imaginary value.
        """
        return self._ctr_data[u'half_width']
    @property
    def quadrature(self):
        """ quadrature (str): Type of quadrature for integrating along the contour.
        """
        return self._ctr_data[u'quad']
    @property
    def energy_min(self):
        """ energy_min (float): minimum (real) energy defining the contour.
        """
        return self._settings[u'energy_min']
    @property
    def energy_max(self):
        """ energy_max (float): maximum (real) energy defining the contour.
        """
        return self._settings[u'energy_max']

    #--------------------------------------------------------------
    def _generateCtrData(self):
        """ Wrapper to populate contour data for a given contour type.
        """
        if   self.name == u'CIRCLE':
            self._ctr_data = self._contourCircle()
        elif self.name == u'CONSTR':
            self._ctr_data = self._contourConstR()
        elif self.name == u'CONSTL':
            self._ctr_data = self._contourConstL()
        elif self.name == u'FERMIS':
            self._ctr_data = self._contourFermiS()
        elif self.name == u'FERMIA':
            self._ctr_data = self._contourFermiA()
        elif self.name == u'EXP':
            self._ctr_data = self._contourExp()
        elif self.name == u'MONOMIAL':
            self._ctr_data = self._contourMonomial()

    #--------------------------------------------------------------
    def updateSettings(self, override):
        """ Update settings with input dictionary and regenerate contour.

        Args:
            override (dict): Settings to override.
        """
        # Override the settings
        for h in override:
            if h not in list(self._settings.keys()):
                raise KeyError(u"Invalid override setting {:} for contour {:}.".format(h,self.name))
            else:
                self._settings[h] = override[h]
        self._generateCtrData()

    #--------------------------------------------------------------
    def setHfbInterval(self, hfb, beta, shift=0.0):
        """
        Update energy interval settings using hfb object properties. This is used
        by default in the pynfma workflow.

        Note:
            Default behavior is to always compute the interval [0, EQRPA_max].
            If (Egs - buffer)<0, compute the larger interval [Egs-buffer,EQRPA_max].
            This will occur if buffer is large, or if Egs<0 for odd nuclei, or both.

        Args:
            hfb (hfbthoRun, dict): The hfb object from which the properties are taken,
                or a dict of the necessary properties.
            beta (str): The beta decay type.
        """

        # Get the relevant HFB quantities from input dict or hfb_obj
        if isinstance(hfb, dict):
            Egs, eqrpamax = hfb[u'E_gs'], hfb[u'EQRPA_max']
            ft_active, temp = hfb[u'ft_active'], hfb[u'temperature']
        else:
            if hfb.beta != beta:
                print(u"WARNING: HFB object beta type missing or different than supplied value.")
                print(u"         Setting attribute 'beta' to supplied value.")
                hfb.beta = beta
            Egs = hfb.soln_dict[u'E_gs']
            eqrpamax = hfb.soln_dict[u'EQRPA_max']
            ft_active = hfb.ft_active
            temp = hfb.temperature

        # Even and odd zero-temp defaults
        buff = self._settings[u'hfb_emin_buff']
        emin = min(Egs-buff, 0.0)
        emax = eqrpamax

        # Finite-temp defaults
        if ft_active:
            # Strength at negative energies is important at FT, but the FT prefactor
            # dies from -inf to 0 at more negative energies. We don't want to include portions
            # where expontential in prefactor overflows and prefactor is exactly zero because
            # the contour integration suffers. The default here is -30, then we override with
            # user input, if any, then we make the contour smaller if necessary in ft_contours().
            emin = -30.0
            # For spontaneous decays we still go to eqrpamax, but for FTEC energetic electrons
            # can be captured and supply the energy to decay even if Q<0, thus emax=+inf here.
            # Here, phase space dies off at large energies, and again will go to zero. Default
            # here is +30, then we override with user input, if any, then make the contour smaller
            # if necessary in ftec_multiple_densities()
            if beta == u'c': emax = 30.0

        self._settings[u'energy_min']= emin + shift
        self._settings[u'energy_max']= emax + shift

        # If nans, avoid errors generating contours, and treat as
        # invalid down the line.
        if np.isnan(Egs) or np.isnan(eqrpamax):
            self._settings[u'energy_min']= 0.0
            self._settings[u'energy_max']= 0.0

        self._generateCtrData()

    #--------------------------------------------------------------
    def _contourCircle(self):
        """
        Generate contour data for:
        Circular contour - centered on real axis, specify left bound, right bound,
        initial theta, and whether to offset from the real axis by shifting the center
        into the complex plane. Equally spaced or gauss-Legendre spaced points available.

        Returns:
            dict
        """

        # Settings
        ctr_bound_left  = self._settings[u'energy_min']
        ctr_bound_right = self._settings[u'energy_max']
        npts       = self._settings[u'nr_points']
        quad       = self._settings[u'beta_quadrature']
        use_gl_ctr = self._settings[u'use_gauleg_ctr']
        shift_imag = self._settings[u'shift_imag']
        t0         = self._settings[u'theta_init']
        max_height = self._settings[u'max_height']

        # EQRPA points
        if use_gl_ctr:
            theta, glwts = np.polynomial.legendre.leggauss(npts)
            t1 = t0+2.0*np.pi
            f1 = 0.5*(t0+t1)
            f2 = 0.5*(t1-t0)
            glwts = glwts*f2
            theta = theta*f2 + f1
        else:
            theta = np.linspace(t0,t0+2.0*np.pi,npts)
            glwts = np.zeros(npts)

        r0 = 0.5*(ctr_bound_right + ctr_bound_left)
        r = ctr_bound_right - r0

        if r > max_height:
            cos = np.cos(theta)
            sin = np.sin(theta)
            ctr_z = r0 + r*cos + max_height*sin*1j
            ctr_dzdt = -r*sin + max_height*cos*1j
        else:
            ctr_z = r0 + r*np.exp(1j*theta)
            ctr_dzdt = 1j*(ctr_z - r0)

        # WE ONLY NEED TO COMPUTE HALF AS MANY POINTS
        nr_compute = ((npts + 1)//2)

        # Shift off real axis. Compute all points
        if shift_imag != 0:
            ctr_z += shift_imag*1j
            nr_compute = npts
        # Rotated contour - Gauleg is only symmetric for t0=integer*pi
        elif use_gl_ctr:
            if abs(np.mod(t0, np.pi)) > 1e-10:
                nr_compute = npts
        # Rotated contour - Equal spacing is only symmetric for t0 a multiple of spacing/2
        else:
            half_spacing = np.pi/(npts-1)
            if abs(np.mod(t0, half_spacing)) > 1e-10:
                nr_compute = npts

        ctr_data = {
                u'nr_points' : npts,   u'nr_compute': nr_compute, u'use_gl_ctr': use_gl_ctr,
                u'ctr_z'     : ctr_z,  u'ctr_dzdt'  : ctr_dzdt,   u'theta'     : theta,
                u'glwts'     : glwts,  u'half_width': None,       u'quad'      : quad,
                u'closed'    : True
                }

        return ctr_data

    #--------------------------------------------------------------
    def _contourConstR(self):
        """
        Generate contour data for:
        A constant line at hw defined via range(min,max,step), where step = hw*de_hw_ratio.
        Spacing is set, and npts is determined.

        Returns:
            dict
        """
        # Settings
        emin = self._settings[u'energy_min']
        emax = self._settings[u'energy_max']
        hw   = self._settings[u'half_width']
        de_hw_ratio = self._settings[u'de_hw_ratio']
        quad = self._settings[u'beta_quadrature']

        # EQRPA points (spacing determined by interval and nr_points)
        de = hw.real*de_hw_ratio
        ctr_z = np.arange(emin, emax+de, de)
        if not np.isreal(hw):
            ctr_z = ctr_z*1j + hw.imag
        else:
            ctr_z = ctr_z + hw.real*1j
        npts = len(ctr_z)

        ctr_data = {
                u'nr_points':npts,        u'nr_compute': npts,        u'use_gl_ctr': False,
                u'ctr_z': ctr_z,          u'ctr_dzdt': np.ones(npts), u'theta': np.zeros(npts),
                u'glwts': np.zeros(npts), u'half_width' : hw,         u'quad' : quad,
                u'closed': False
                }

        return ctr_data

    #--------------------------------------------------------------
    def _contourConstL(self):
        """
        Generate contour data for:
        A constant line at hw, define via linspace(min,max,npts).
        npts is set, spacing is determined.

        Returns:
            dict
        """
        # Settings
        emin = self._settings[u'energy_min']
        emax = self._settings[u'energy_max']
        npts = self._settings[u'nr_points']
        hw   = self._settings[u'half_width']
        quad = self._settings[u'beta_quadrature']

        # EQRPA points (spacing determined by interval and nr_points)
        ctr_z = np.linspace(emin, emax, npts)
        if not np.isreal(hw):
            ctr_z = ctr_z*1j + hw.imag
        else:
            ctr_z = ctr_z + hw.real*1j

        ctr_data = {
                u'nr_points':npts,        u'nr_compute': npts,        u'use_gl_ctr': False,
                u'ctr_z': ctr_z,          u'ctr_dzdt': np.ones(npts), u'theta': np.zeros(npts),
                u'glwts': np.zeros(npts), u'half_width' : hw,         u'quad' : quad,
                u'closed': False
                }

        return ctr_data

    #--------------------------------------------------------------
    def _contourFermiA(self):
        """
        Generate contour data for:
        Adaptive Fermi Fct. npts is set, spacing is determined, except in extreme cases.

        The adaptive fermi contour changes the contour itself to fill an energy interval
        with a specified number of points, prioritizing small width at small energy.
        This occurs in stages, between two limiting cases.

        1. Max points limiting case: CONSTL at a min half-width.
        2. Raise right side of contour to form a fermi function type profile.
        3. If right side reaches max half-width, raise left side of contour until max half-width.
        4. Min points limiting case: CONSTL at max half-width
        5. If the interval is still too large after stage 4, the number of points is increased until
           CONSTL at max half-width fills the interval.

        Returns:
            dict
        """
        # Settings
        npts_min = self._settings[u'nr_points_min']
        npts_cap = self._settings[u'nr_points_max']
        emin_in  = self._settings[u'energy_min']
        emax_in  = self._settings[u'energy_max']
        hwmin    = self._settings[u'hw_min']
        hwmax    = self._settings[u'hw_max']
        u_i      = self._settings[u'u_percent_interval']
        t_i      = self._settings[u't_percent_interval']
        de_hw_ratio = abs(self._settings[u'de_hw_ratio'])
        quad     = self._settings[u'beta_quadrature']

        shift_increment = 0.0001

        # Convert interval to start at 0 so we take the appropriate part of the fermi fct
        emin = float(emin_in - emin_in)
        emax = float(emax_in - emin_in)
        de_min = float(hwmin*de_hw_ratio)
        de_max = float(hwmax*de_hw_ratio)

        # If the interval is too large, we are forced to increase npts (note wmax = interval)
        abs_npts_min = int(np.ceil( (emax/de_max) )) + 1
        if npts_cap <= abs_npts_min:
            print(u"WARNING: Absolute min npts is {:}, but requested {:}.".format(abs_npts_min, npts_cap))
            print(u"         Increasing npts to fill line at max half width.")
            w = np.arange(emin, emax+de_max, de_max)
            ctr_z = (w + emin_in) + (np.ones(len(w))*hwmax)*1j
        else:
            # The Fermi function
            u = u_i*emax
            t = t_i*emax
            def fermi(x):
                return 1.0 - 1.0/(np.exp(((x - u)/t)) + 1.0)

            # Construct the contour
            shift_u, shift_l = 0.0, 0.0
            shift_max = hwmax - hwmin
            while True:
                # Get the energy grid knowing function for de
                if shift_u == 0.0 and shift_l == 0.0:
                    fct = 0
                    w = np.arange(emin, emax+de_min, de_min)
                    if len(w) < npts_min: w = np.linspace(emin, emax, npts_min)
                elif shift_u >= shift_max and shift_l >= shift_max:
                    fct = 2
                    w = np.arange(emin, emax+de_max, de_max)
                else:
                    fct = 1
                    w = np.zeros(npts_cap); w[0] = emin
                    for i in range(npts_cap-1):
                        # negative de_hw_ratio will set de based on arc length
                        # turns out since half width is so small compared to de this doesn't
                        # matter much (and takes much longer to converge for some reason)
                        if self._settings[u'de_hw_ratio'] < 0:
                            darc= (a*fermi(w[i])+b)*de_hw_ratio
                            dhw = (a*fermi(w[i])+b) - (a*fermi(w[i-1])+b)
                            de_calc  = np.sqrt(darc**2 - dhw**2)
                        else:
                            de_calc = (a*fermi(w[i])+b)*de_hw_ratio
                        w[i+1] = w[i] + de_calc


                # Adjust function until we satisfy constraints
                if w[-1] >= emax and len(w) <= npts_cap:
                    break
                else:
                    if shift_u < shift_max:
                        shift_u += shift_increment
                    else:
                        shift_u = shift_max
                        shift_l += shift_increment
                        if shift_l > shift_max:
                            shift_l = shift_max
                    ymax = hwmin + shift_u
                    ymin = hwmin + shift_l
                    a = ((ymax-ymin)/(fermi(emax)- fermi(emin)))
                    b = ymin - a*fermi(emax)

            # Fudge the endpoint to include emax
            if w[-1] != emax: w[-1] = emax

            # Return the contour on the energy grid
            if fct == 0:
                ctr_z = (w + emin_in) + (np.ones(len(w))*hwmin)*1j
            elif fct == 2:
                ctr_z = (w + emin_in) + (np.ones(len(w))*hwmax)*1j
            else:
                ctr_z = (w + emin_in) + (a*fermi(w) + b)*1j

        npts = len(ctr_z)

        ctr_data = {
                u'nr_points': npts,       u'nr_compute': npts,        u'use_gl_ctr': False,
                u'ctr_z': ctr_z,          u'ctr_dzdt': np.ones(npts), u'theta': np.zeros(npts),
                u'glwts': np.zeros(npts), u'half_width' : None,       u'quad' : quad,
                u'closed':False
                }

        return ctr_data

    #--------------------------------------------------------------
    def _contourFermiS(self):
        """
        Generate contour data for:
        Static fermi fct -- npts are determined not set.

        Returns:
            dict
        """
        # Settings
        emin_in  = self._settings[u'energy_min']
        emax_in  = self._settings[u'energy_max']
        hwmin    = self._settings[u'hw_min']
        hwmax    = self._settings[u'hw_max']
        u_i      = self._settings[u'u_percent_interval']
        t_i      = self._settings[u't_percent_interval']
        de_hw_ratio = abs(self._settings[u'de_hw_ratio'])
        quad     = self._settings[u'beta_quadrature']

        # Convert interval to start at 0 so we take the appropriate part of the fermi fct
        emin = float(emin_in - emin_in)
        emax = float(emax_in - emin_in)
        de_min = float(hwmin*de_hw_ratio)

        abs_npts_max = int(np.ceil((emax/de_min))) + 1

        # The Fermi function
        u = u_i*emax
        t = t_i*emax
        def func(x):
            return 1.0 - 1.0/(np.exp(((x - u)/t)) + 1.0)
        a = ((hwmax-hwmin)/(func(emax) - func(emin)))
        b = hwmax - a*func(emax)

        # Get the energy grid knowing function for de
        w    = np.zeros(abs_npts_max)
        w[0] = emin
        i=0
        while True:
            de_calc = (a*func(w[i])+b)*de_hw_ratio
            w[i+1] = w[i] + de_calc
            if w[i+1] >= emax:
                break
            else:
                i += 1
        w = np.trim_zeros(w, u'b')

        # Fudge the endpoint to include emax
        w[-1] = emax

        # Return the contour on the energy grid
        ctr_z = (w + emin_in) + (a*func(w) + b)*1j

        npts = len(ctr_z)

        ctr_data = {
                u'nr_points': npts,       u'nr_compute': npts,        u'use_gl_ctr': False,
                u'ctr_z': ctr_z,          u'ctr_dzdt': np.ones(npts), u'theta': np.zeros(npts),
                u'glwts': np.zeros(npts), u'half_width' : None,       u'quad' : quad,
                u'closed':False
                }

        return ctr_data

    #--------------------------------------------------------------
    def _contourExp(self):
        """
        Generate contour data for:
        Static exponential -- npts are determined not set.

        Returns:
            dict
        """
        # Settings
        emin_in  = self._settings[u'energy_min']
        emax_in  = self._settings[u'energy_max']
        hwmin    = self._settings[u'hw_min']
        hwmax    = self._settings[u'hw_max']
        p_i      = self._settings[u'p_percent_interval']
        de_hw_ratio = abs(self._settings[u'de_hw_ratio'])
        quad     = self._settings[u'beta_quadrature']

        # Convert interval to start at 0 so we take the appropriate part of the fct
        emin = float(emin_in - emin_in)
        emax = float(emax_in - emin_in)
        de_min = float(hwmin*de_hw_ratio)
        abs_npts_max = int(np.ceil((emax/de_min))) + 1

        # Define the fct
        p = p_i*emax
        def func(x):
            return np.exp((x/p))
        a = ((hwmax - hwmin)/(func(emax) - func(emin)))
        b = hwmin - a*func(emin)

        # Calculate x values to fill the interval
        w    = np.zeros(abs_npts_max)
        w[0] = emin; i=0
        while True:
            de_calc = (a*func(w[i])+b)*de_hw_ratio
            w[i+1] = w[i] + de_calc
            if w[i+1] >= emax:
                break
            else:
                i += 1
        w = np.trim_zeros(w, u'b')

        # Fudge the endpoint to include emax
        w[-1] = emax

        # Return the contour on the energy grid
        ctr_z = (w + emin_in) + (a*func(w) + b)*1j

        npts = len(ctr_z)

        ctr_data = {
                u'nr_points': npts,       u'nr_compute': npts,        u'use_gl_ctr': False,
                u'ctr_z': ctr_z,          u'ctr_dzdt': np.ones(npts), u'theta': np.zeros(npts),
                u'glwts': np.zeros(npts), u'half_width' : None,       u'quad' : quad,
                u'closed':False
                }

        return ctr_data

    #--------------------------------------------------------------
    def _contourMonomial(self):
        """
        Generate contour data for:
        Static Monomial (power>0) -- npts are determined not set.

        Returns:
            dict
        """
        # Settings
        emin_in  = self._settings[u'energy_min']
        emax_in  = self._settings[u'energy_max']
        hwmin    = self._settings[u'hw_min']
        hwmax    = self._settings[u'hw_max']
        power    = self._settings[u'power']
        if power <= 0.0:
            raise ValueError(u"Monomial power must be >= 0")
        de_hw_ratio = abs(self._settings[u'de_hw_ratio'])
        quad     = self._settings[u'beta_quadrature']

        # Convert interval to start at 0 so we take the appropriate part of the fct
        emin = float(emin_in - emin_in)
        emax = float(emax_in - emin_in)
        de_min = float(hwmin*de_hw_ratio)
        abs_npts_max = int(np.ceil((emax/de_min))) + 1

        # Define the fct
        def func(x):
            return x**power
        a = ((hwmax - hwmin)/(func(emax) - func(emin)))
        b = hwmin - a*func(emin)

        # Calculate x values to fill the interval
        w    = np.zeros(abs_npts_max)
        w[0] = emin; i=0
        while True:
            de_calc = (a*func(w[i])+b)*de_hw_ratio
            w[i+1] = w[i] + de_calc
            if w[i+1] >= emax:
                break
            else:
                i += 1
        w = np.trim_zeros(w, u'b')

        # Fudge the endpoint to include emax
        w[-1] = emax

        # Return the contour on the energy grid
        ctr_z = (w + emin_in) + (a*func(w) + b)*1j

        npts = len(ctr_z)

        ctr_data = {
                u'nr_points': npts,       u'nr_compute': npts,        u'use_gl_ctr': False,
                u'ctr_z': ctr_z,          u'ctr_dzdt': np.ones(npts), u'theta': np.zeros(npts),
                u'glwts': np.zeros(npts), u'half_width' : None,       u'quad' : quad,
                u'closed':False
                }

        return ctr_data

    #-----------------------------------------------------------------------
    def _getDefaults(self):
        """
        Return the dictionary for default contour settings.

        Returns:
            dict

        Raises:
            ValueError
        """

        default_raw = copy.deepcopy(DEFAULTS[u'ctr'])

        # Get the common "interval" settings, plus those for the specific contour
        default_settings = default_raw[u'INTERVAL']
        try:
            default_settings.update(default_raw[self.name])
        except KeyError:
            raise ValueError(u"Requested contour {:} not implemented.".format(self.name))

        return default_settings
