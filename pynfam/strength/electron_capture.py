# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
# -------------- Utilities -----------------
import numpy as np
from scipy.optimize import brentq
from scipy.special import gamma
# ----------- Relative Imports -------------
from ..config import HBAR_MEC, MEC2, KB, NA, ALPHA, R0
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2020-11-06'


#=======================================================================
#             Terrestrial electron capture from shell 'X'
#=======================================================================
def kappa_x(n,l,j):
    """ Dirac angular momentum quantum number """
    if j == 2*l - 1: return l
    if j == 2*l + 1: return -(l+1)

#-----------------------------------------------------------------------
def betasq_x(n,l,j,Z,A):
    """ Hydrogen-like wavefunction evaluated at nuclear surface (squared).
    Also known as the coulomb amplitude """
    if (n,l,j) != (1,0,1):
        raise ValueError("Only K-shell coulomb amplitudes are implemented")

    if (n,l,j) == (1,0,1):
        Bk = energy_x(n,l,j,Z)
        Zk = Z - 0.35 # Effective K-shell charge due to screening of other electrons
        s = 1 - Bk
        R = R0/HBAR_MEC*A**(1.0/3.0)
        g_1k_2 = (2-Bk)/(2-gamma(2*s+1))*(2*ALPHA*Zk)**3*(2*ALPHA*Zk*R)**(2*s-2)
        return g_1k_2

#-----------------------------------------------------------------------
def energy_x(n,l,j,Z):
    """ Hydrogen-like electron binding energy """
    if (n,l,j) != (1,0,1):
        raise ValueError("Only K-shell binding energies are implemented")

    if (n,l,j) == (1,0,1):
        Zk = Z - 0.35 # Effective K-shell charge due to screening of other electrons
        return 1 - np.sqrt(1 - (ALPHA*Zk)**2)


#=======================================================================
#               Stellar electron capture at finite temp
#=======================================================================
def gen_fermi_dirac(x, k, eta, beta):
    """ Generalized Fermi-Dirac integrand. """
    def f(x):
        return x**k * np.sqrt(1.0 + 0.5*x*beta) / (np.exp(x-eta) + 1.0)
    def f_over(x):
        return x**k * np.sqrt(1.0 + 0.5*x*beta) * np.exp(eta-x)
    try:
        return np.array([f(xx) if (xx-eta) < 100 else f_over(xx) for xx in x])
    except TypeError:
        return np.array([f(xx) if (xx-eta) < 100 else f_over(xx) for xx in [x]])

#-----------------------------------------------------------------------
def gen_fermi_dirac_int(k, eta, beta, gpts = 20):
    """ Generalized Fermi-Dirac integral.

       https://iopscience.iop.org/article/
       10.1086/313121/fulltext/37118.text.html
       J.M. Aparicio, Astro phys J Supp Series 117:627-632, 1998)

       Checked against mathematica -- Looks great
    """
    # Gauss-legendre shortcut routine
    def transform_interval(x, wts, xmin, xmax):
        """ Transform Gauss-Legendre abcissae and weights for integration from xmin to xmax. """
        f1 = 0.5*(xmin+xmax)
        f2 = 0.5*(xmax-xmin)
        wts = wts*f2
        x = x*f2 + f1
        return x, wts

    # Optimal parameters (obtained by authors via simulated annealing)
    D  = 3.3609; sig= 9.1186e-2
    a1 = 6.7774; b1 = 1.1418   ; c1 = 2.9826
    a2 = 3.7601; b2 = 9.3719e-2; c2 = 2.1063e-2; d2 = 3.1084e1; e2 = 1.0056
    a3 = 7.5669; b3 = 1.1695   ; c3 = 7.5416e-1; d3 = 6.6558  ; e3 = -1.2819e-1

    arg = sig*(eta - D)
    if (arg <= 5e1):
        xi = np.log(1.0 + np.exp(arg))/sig
    else:
        xi = eta - D
    Xa = (a1 + b1*xi +    c1*xi*xi)/(1 + c1*xi)
    Xb = (a2 + b2*xi + c2*d2*xi*xi)/(1 + e2*xi + c2*xi*xi)
    Xc = (a3 + b3*xi + c3*d3*xi*xi)/(1 + e3*xi + c3*xi*xi)

    S1 = Xa - Xb
    S2 = Xa
    S3 = Xa + Xc

    x, gwts = np.polynomial.legendre.leggauss(gpts)

    # 1st interval - Legendre with x = z^2 (x_Legendre = z, glwts ~ dz)
    x1, gwts1 = transform_interval(x, gwts, 0, np.sqrt(S1))
    gwts1 = 2.0*x1*gwts1 # include 2*z*dz
    f1 = np.sum(gen_fermi_dirac(x1*x1, k, eta, beta)*gwts1)

    # 2nd interval - Legendre
    x2, gwts2 = transform_interval(x, gwts, S1, S2)
    f2 = np.sum(gen_fermi_dirac(x2, k, eta, beta)*gwts2)

    # 3rd interval - Legendre
    x3, gwts3 = transform_interval(x, gwts, S2, S3)
    f3 = np.sum(gen_fermi_dirac(x3, k, eta, beta)*gwts3)

    # 4th interval - Laguerre with x = z + S3
    x, gwts = np.polynomial.laguerre.laggauss(gpts)
    x4    = x + S3
    gwts4 = gwts*np.exp(x) # Include inverse of glag weight
    f4 = np.sum(gen_fermi_dirac(x4, k, eta, beta)*gwts4)

    return f1+f2+f3+f4

#-----------------------------------------------------------------------
def charge_neutrality(uf, Ye_rho, T):
    """ Charge neutrality condition for electron-positron Fermi gas. """

    pre = np.sqrt(2)/((HBAR_MEC*1e-13)**3*np.pi**2*NA) # dimensionful prefactor
    eta = uf/(KB*T)
    beta= (KB*T)/MEC2
    etap= - eta - 2.0/beta

    FD_ele = gen_fermi_dirac_int(0.5, eta,  beta) + beta*gen_fermi_dirac_int(1.5, eta,  beta)
    if beta > 0.02:
        FD_pos = gen_fermi_dirac_int(0.5, etap, beta) + beta*gen_fermi_dirac_int(1.5, etap, beta)
    else:
        FD_pos = 0.0

    return Ye_rho - pre*beta**(1.5)*(FD_ele - FD_pos)

#-----------------------------------------------------------------------
def chem_pot(logYe_rho, T):
    """ Solve charge neutrality equation using root finder to determine electron chemical potential.
    Returns chemical potential in MeV """

    # Allow iteration over scalars
    logYe_rho = np.array(logYe_rho, copy=False, ndmin=1)

    TK = T/KB # Temp in Kelvin
    Ye_rho = 10**logYe_rho # undo the log

    # chemical potential approaches -mec2 for high temps. Largest values occur for
    # high densities and low temps, no more than 1e4...
    return np.array([brentq(lambda uf : charge_neutrality(uf, Ye_rho=r, T=TK), -.6, 1e10) for r in Ye_rho])

#-----------------------------------------------------------------------
def chem_pot_0T(logYe_rho):
    """ Analytic expression for chemical potential at zero temperature for any density.
    Returns chemical potential in MeV
    """

    return MEC2*(np.sqrt((1.02e-4)*(10**logYe_rho)**(2/3) + 1) - 1)
