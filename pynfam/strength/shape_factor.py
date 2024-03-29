# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import range
from builtins   import object
from builtins   import str
# -------------- Utilities -----------------
from collections     import defaultdict
from scipy.integrate import trapz, simps, cumtrapz
import numpy  as np
import pandas as pd
import datetime
import os
# ----------- Relative Imports -------------
from ..config        import HBAR_MEC, MEC2, MN, KAPPA, ALPHA, R0
from .fam_strength   import famStrength
from .phase_space    import phaseSpace, mu_ke, gam_ke
from ..utilities.hfb_utils import convString
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#=========================================================================#
#                         Class shapeFactor                               #
#=========================================================================#
class shapeFactor(phaseSpace):
    """
    A shapeFactor is a phaseSpace subclass that contains a set of famStrength objects
    and their corresponding contour, and the methods to compute the phase-space-
    weighted shape factors from the raw strengths for allowed and first-forbidden
    beta decay contributions. It also contains the methods to integrate the shape
    factor for beta decay rates, and write the results to file.

    Args:
        strengths (list of famStrength): The strengths to be combined.

    Attributes:
        strengths (dict): Contains {s.genopname:s for s in strengths}
        contour (famContour): The contour object for the strengths.
        nucleus (tuple): (N,Z,A)
        ops (list of str): A list of genopnames for the strengths.
        beta (str): The beta decay type.
        is_raw (bool): Indicates whether the shape factor was computed
            with (True) or without (False) adjusted strength functions.
            (e.g. with cstr_df_zeroed).
        sf_df (dataframe): The shape factor along the contour, split by contribution.
        ps_df (dataframe): The different phase space integrals along the contour.
        ps_energy_df (dataframe): The EQRPA and W0 for the phase space integrals.
        sf_metadict (dict): Various values to print to the output header.
        e_1stpeak (float): Energy of the first peak in the strength indicating the
            QRPA estimate of Egs. None for closed contours.
        ind_1stpeak (int): Index of the location of e_1stpeak. None for closed contours.
        ind_start (int): Index of the location of the base of the first peak.
            None for closed contours.
    """

    # Everything we need to integrate         # open ctr      closed ctr
    #sf_df[u'EQRPA'] = eqrpa                  # real(EQRPA)   complex(EQRPA)
    #sf_df[u'dzdt']  = self.contour.ctr_dzdt  # np.ones(dim)  complex(dz/dt)
    #sf_df[u'theta'] = self.contour.theta     # np.zeros(dim) float(theta)
    #sf_df[u'glwts'] = self.contour.glwts     # np.zeros(dim) float(glwts)

    # Complete list of keys for rate contributions in output files
    beta_totals      = [u'Total'                ,
                        u'Total-Allowed'        ,
                        u'Total-GT'             ,
                        u'Total-Forbidden'
                        ]
    beta_al_contribs = [u'Allowed-Fermi'        ,
                        u'Allowed-GT_K=0'       ,
                        u'Allowed-GT_K=1'
                        ]
    beta_ffK_totals  = [u'Forbidden-K=0'        ,
                        u'Forbidden-K=1'        ,
                        u'Forbidden-K=2'
                        ]
    beta_ffJ_totals  = [u'Forbidden-J=0'        ,
                        u'Forbidden-J=1'        ,
                        u'Forbidden-J=2'
                        ]
    beta_ff_contribs = [u'Forbidden-(J,K)=(0,0)',
                        u'Forbidden-(J,K)=(1,0)',
                        u'Forbidden-(J,K)=(1,1)',
                        u'Forbidden-(J,K)=(2,0)',
                        u'Forbidden-(J,K)=(2,1)',
                        u'Forbidden-(J,K)=(2,2)'
                        ]
    betaout_keys = beta_totals + beta_al_contribs + beta_ffK_totals + beta_ffJ_totals + beta_ff_contribs


    """ Labels for the different beta decay contributions used in output files.
    """

    # Complete list of 14 operators by K (genop to ignore beta type)
    al_ops_raw  = [u'F_K0', u'GT_K0', u'GT_K1']
    k0_fops_raw = [u'P_K0', u'R_K0', u'PS0_K0', u'RS0_K0', u'RS1_K0', u'RS2_K0']
    k1_fops_raw = [u'P_K1', u'R_K1', u'RS1_K1', u'RS2_K1']
    k2_fops_raw = [u'RS2_K2']
    #11/26/23: list of operators by J + K 
    j0_k0_fops_raw = [u'PS0_K0', u'RS0_K0']
    j1_k0_fops_raw = [u'P_K0', u'R_K0', u'RS1_K0']
    j1_k1_fops_raw = [u'P_K1', u'R_K1', u'RS1_K1']
    fops_j_dict = {u'PS0': 0, u'RS0': 0, u'P': 1, u'R': 1, u'RS1': 1, u'RS2': 2}
    # Complete list of 8 bare operators and 4 cross terms
    al_ops = [u'F', u'GT']
    ff_ops = [u'P', u'R', u'PS0', u'RS0', u'RS1', u'RS2']
    xts    = [u'RS0_PS0', u'R_RS1', u'P_RS1', u'R_P']

    # Shape factor labels by J. The R index indicates the phase space factor to use
    CJ0_keys = [u'J0_R1',u'J0_R2']
    CJ1_keys = [u'J1_R1',u'J1_R2',u'J1_R3',u'J1_R4',u'J1_R5',u'J1_R6']
    CJ2_keys = [u'J2_R5',u'J2_R6']
    C_keys   = CJ0_keys+CJ1_keys+CJ2_keys

    #--------------------------------------------------------------
    def __init__(self, strengths, beta_type, ps_contour=None):
        if not isinstance(strengths, list):
            strengths = [strengths]
        if not all([isinstance(s, famStrength) for s in strengths]):
            raise ValueError(u"Strength inputs must be famStrength objects.")

        phaseSpace.__init__(self, beta_type)

        # Some consistency checks for the provided strengths
        ops = [s.genopname for s in strengths]
        nuc = strengths[0].nucleus
        dim = strengths[0].contour.nr_points
        sbeta= strengths[0].beta
        if any([s.nucleus != nuc for s in strengths]): raise ValueError
        if any([s.str_df.shape[0] != dim for s in strengths]): raise ValueError
        if any([len(s.contour.theta)    != dim for s in strengths]): raise ValueError
        if any([len(s.contour.glwts)    != dim for s in strengths]): raise ValueError
        if any([len(s.contour.ctr_z)    != dim for s in strengths]): raise ValueError
        if any([len(s.contour.ctr_dzdt) != dim for s in strengths]): raise ValueError
        if any([s.beta != sbeta for s in strengths]): raise ValueError

        # Consistency check on beta type between strength and phase space
        if beta_type in [u'+',u'c']:
            if sbeta != u'+': raise ValueError
        else:
            if sbeta != u'-': raise ValueError

        # Set some main attributes
        self.strengths = {s.genopname : s for s in strengths}
        self.contour   = strengths[0].contour
        self.nucleus   = nuc
        self.ops       = ops
        self.is_raw    = False
        self.sf_metadict= {u'beta_type':self.beta, u'quadratr':None, u'FAM_ctr':self.contour.name_and_int}
        if not ps_contour:
            self.ps_contour = self.contour
        else:
            if ps_contour.closed != self.contour.closed:
                raise RuntimeError("Alternate phase space contour type (open/closed) should match.")
            self.ps_contour = ps_contour

        # Phase space weighted shape factor solution
        self.sf_df       = None
        self.ps_df       = None
        self.ps_energy_df= None
        self.e_1stpeak   = None
        self.ind_1stpeak = None
        self.ind_start   = None
        self._ps_test    = None

    # Split real and imaginary parts because pandas doesn't do well with complex
    @property
    def sf_df_real(self):
        """ sf_df_real (dataframe): The real part of the shape factor contributions.
        """
        if self.sf_df is None: return None
        re_sf = pd.DataFrame()
        re_sf[u'Re(EQRPA)'] = np.real(self.contour.ctr_z)
        re_sf[u'Im(EQRPA)'] = np.imag(self.contour.ctr_z)
        for col in self.sf_df:
            re_sf[col] = np.real(self.sf_df[col].values)
        return re_sf

    @property
    def sf_df_imag(self):
        """ sf_df_imag (dataframe): the imaginary part of the shape factor contributions.
        """
        if self.sf_df is None: return None
        im_sf = pd.DataFrame()
        im_sf[u'Re(EQRPA)'] = np.real(self.contour.ctr_z)
        im_sf[u'Im(EQRPA)'] = np.imag(self.contour.ctr_z)
        for col in self.sf_df:
            im_sf[col] = np.imag(self.sf_df[col].values)
        return im_sf

    @property
    def ps_df_real(self):
        """
        ps_df_real (dataframe): the real part of the phase space contributions.
        W0 is the maximum electron energy for a given transition with energy Re(EQRPA):

        :math:`W_0 = 1 + (EQRPA_{max} - EQRPA)/m_e c^2`
        """
        if self.ps_df is None: return None
        re_ps= pd.DataFrame()
        re_ps[u'Re(EQRPA)'] = self.ps_energy_df[u'Re(EQRPA)']
        re_ps[u'Im(EQRPA)'] = self.ps_energy_df[u'Im(EQRPA)']
        re_ps[u'Re(W0)']    = self.ps_energy_df[u'Re(W0)']
        re_ps[u'Im(W0)']    = self.ps_energy_df[u'Im(W0)']
        for col in self.ps_df:
            re_ps[col] = np.real(self.ps_df[col].values)
        return re_ps

    @property
    def ps_df_imag(self):
        """
        ps_df_imag (dataframe): the imaginary part of the phase space contributions.
        W0 is the maximum electron energy for a given transition with energy Re(EQRPA):

        :math:`W_0 = 1 + (EQRPA_{max} - EQRPA)/m_e c^2`
        """
        if self.ps_df is None: return None
        im_ps = pd.DataFrame()
        im_ps[u'Re(EQRPA)'] = self.ps_energy_df[u'Re(EQRPA)']
        im_ps[u'Im(EQRPA)'] = self.ps_energy_df[u'Im(EQRPA)']
        im_ps[u'Re(W0)']    = self.ps_energy_df[u'Re(W0)']
        im_ps[u'Im(W0)']    = self.ps_energy_df[u'Im(W0)']
        for col in self.ps_df:
            im_ps[col] = np.imag(self.ps_df[col].values)
        return im_ps

    @property
    def genops(self):
        """ genops (list of str): The keys of the strengths attribute, Bare_operator+_+K.
        """
        return list(self.strengths.keys())

    @property
    def sf_metastr(self):
        """ sf_metastr (str): The header string for output files, constructed from sf_metadict.
        """
        if self.sf_metadict is None: return None
        return self.sfSummaryString(self.sf_metadict)

    @property
    def ffk_err(self):
        """ ffk_err (list): If any FF operators missing for a given K=0,1,2, contains True.

            Check we have all forbidden operators for a given K to avoid computing
            weird rates. This could be improved by separating everything by J since
            we do not have cross terms between operators with different J.
        """
        ffk_err = [any([o not in self.genops for o in fops]) \
            for fops in [shapeFactor.k0_fops_raw, shapeFactor.k1_fops_raw, shapeFactor.k2_fops_raw]]
        return ffk_err
    
    #11/26/23: condition for 
    def ffjk_err(self, op, K):
        """ ffjk_err (bool): if any FF operators with crossterms are missing for a given J+K, contains True.
            
            ex: if PS0_K0 is requested but RS0_K0 is not in genops, then return True (no calculations). 
            If P_0 is requested, then R_0 and RS1_0 must be in genops as well. 
            Since RS2 has no crossterms, this will always return False for that operator.
        """
        j_value = shapeFactor.fops_j_dict[op]
        if j_value == 0:
            ffjk_err = any([o not in self.genops for o in shapeFactor.j0_k0_fops_raw])
        elif j_value == 1:
            if K == 0:
                ffjk_err = any([o not in self.genops for o in shapeFactor.j1_k0_fops_raw])
            else:
                ffjk_err = any([o not in self.genops for o in shapeFactor.j1_k1_fops_raw])
        else:
            return False #if J=2, don't need to worry about cross-terms.
        return ffjk_err

    @property
    def all_conv(self):
        """ all_conv (str): Indicates whether all points in all strengths are converged
        """
        try:
            conv_list = [s.meta[u'Conv'] for s in list(self.strengths.values())]
            return convString(conv_list)
        except:
            return u"Error"

    @property
    def total_time(self):
        """ total_time (float): Total run time for all points in all strengths
        """
        try:
            return sum([s.meta[u'Time'] for s in list(self.strengths.values())])
        except:
            return np.nan

    #---------------------------------------------------------------------------
    def _calcShapeFactor(self, hfb_gs, zero_neg, Evloss=False):
        """
        Execute all the methods to combine strengths and phase space into the
        shape factor. See (Mustonen et al., Phys. Rev. C 90, 024308 (2014)).

        Args:
            hfb_gs (hbfthoRun): The hfbthoRun object for the ground state.
            zero_neg (bool): Used the adjusted strength with zeroed negative components.
        """

        self.prepConstants(hfb_gs)
        self.calcPhaseSpace(Evloss)
        b = self.prepStrengths(zero_neg)
        c = self.calcSfByJ(b)
        self.calcSfContributions(c)

        # Get location and index of 1st peak. Default is None for closed ctrs
        if not self.contour.closed and not zero_neg:
            total_str = np.imag(self.sf_df[u'Total'].values)
            eqrpa     = np.real(self.contour.ctr_z)
            ind_1stpk, ind_start = self.findFirstPeak(eqrpa, total_str)
            self.ind_start   = ind_start
            self.ind_1stpeak = ind_1stpk
            self.e_1stpeak   = eqrpa[ind_1stpk]

        # If we found the 1st peak, use this as fam_Qval, otherwise default to hfb
        if self.e_1stpeak is None:
            fam_qval = None
        else:
            fam_qval = self.sf_metadict[u"EQRPAmax"] - self.e_1stpeak

        # Store peak info as meta data
        self.sf_metadict.update({u'E_1stPeak':self.e_1stpeak, u'FAM_Qval':fam_qval})

        # If we've altered the strength, make it known
        if not zero_neg or self.contour.closed:
            self.is_raw = True
        else:
            self.is_raw = False

    #--------------------------------------------------------------
    def calcShapeFactor(self, hfb_gs, zero_neg=False, Evloss=False):
        """
        Public wrapper for _calcShapeFactor using the raw strength as it
        comes out of pnfam. (zero_neg=False).

        Args:
            hfb_gs (hfbthoRun): The ground state hfbthoRun object.
        """
        if zero_neg and self.contour.closed:
            print(u"WARNING: Do not set negative strength to zero for closed contours. Skipping.")
            zero_neg = False

        if not zero_neg:
            self._calcShapeFactor(hfb_gs, zero_neg=False, Evloss=Evloss)
        else:
            # We need the raw strength to estimate where it starts - the 1st positive peak
            if not self.is_raw:
                self._calcShapeFactor(hfb_gs, zero_neg=False, Evloss=Evloss)
            ind_start = self.ind_start

            # Now calculate the zeroed strength (skips looking for 1st peak here)
            self._calcShapeFactor(hfb_gs, zero_neg=True, Evloss=Evloss)
            if ind_start > 0: self.sf_df.loc[:ind_start] = 0.0

    #---------------------------------------------------------------------------
    def prepConstants(self, hfb_gs):
        """ Computes and gathers constants needed to construct shape factor,
        stores them in a dictionary.
        """

        if hfb_gs.beta != self.beta:
            print(u"WARNING: HFB object beta type missing or different than supplied value.")
            print(u"         Setting attribute 'beta' to supplied value.")
            hfb_gs.beta = self.beta

        hfb_qval = hfb_gs.soln_dict[u'HFB_Qval']
        eqrpamax = hfb_gs.soln_dict[u'EQRPA_max']
        LAM = self.LAM
        GA = self.GA
        GV = self.GV
        A = self.nucleus[2]
        rad = (R0*A**(1.0/3.0)/HBAR_MEC)

        Zi = self.nucleus[1]
        if self.beta == u'-':
            Zf = Zi + 1 # Final nucleus
            Zd = Zf     # Z for coulomb distortion - daughter
            LAMd = LAM  # Axial vector coupling
            ec = +1     # Matrix element sign changes
        elif self.beta == u'+':
            Zf = Zi - 1 # Final nucleus
            Zd = -Zf    # Z for coulomb distortion - daughter, -Z for positron wf
            LAMd = -LAM # Axial vector coupling - includes matrix element sign changes
            ec = +1     # Matrix element sign changes
        elif self.beta == u'c':
            # (treat separately from beta+ b/c different energy regimes? Only for closed ctrs?)
            Zf = Zi - 1 # Final nucleus
            Zd = Zi     # Z for coulomb distortion - parent
            LAMd = LAM  # Axial vector coupling
            ec = -1     # Matrix element sign changes

        # Half_width (if it's well defined for the contour)
        half_width = np.imag(self.contour.ctr_z[0])
        if self.contour.closed or any(hw != half_width for hw in np.imag(self.contour.ctr_z)):
            half_width = None

        # Adjust Q if requested (do it here so all values with Q_eff are printed to header)
        Q_eff  = self._settings[u'Q_eff']
        Q_mode = self._settings[u'Q_eff_mode']
        if Q_eff is not None and Q_mode != 0:
            Egs = eqrpamax - hfb_qval
            # Mode 1: q_eff replaces q_hfb
            if abs(Q_mode) == 1:
                hfb_qval = Q_eff
            # Mode 2: q_eff is a shift from q_hfb
            elif abs(Q_mode) == 2:
                hfb_qval = hfb_qval + Q_eff
            eqrpamax = hfb_qval + Egs

        # Store constants
        self.sf_metadict.update({\
                   u'ec'         : ec,   # To relate SF b- to bound ec
                   u'Z'          : Zd,   # To relate SF b- to b+
                   u'LAM'        : LAMd, # To relate SF b- to b+
                   u'ft_active'  : hfb_gs.ft_active,
                   u'temper'     : hfb_gs.temperature,
                   u'HFB_Qval'   : hfb_qval,
                   u'EQRPAmax'   : eqrpamax,
                   u'Half_Width' : half_width,
                   u'A'          : A,
                   u'Zi'         : Zi,
                   u'Zf'         : Zf,
                   u'|gA|/gV'    : LAM,
                   u'gA'         : GA,
                   u'gV'         : GV,
                   u'M_nucleon'  : MN,
                   u'Radius'     : rad,
                   u'alpha*Z'    : ALPHA*abs(Zd),
                   u'alpha*Z/2R' : (ALPHA*abs(Zd)/(2.0*rad)),
                   u'W0_max'     : 1.0 + (hfb_qval/MEC2),
                   u'W0*R'       : ((1.0+hfb_qval)/MEC2*rad)
                   })

    #---------------------------------------------------------------------------
    def calcPhaseSpace(self, Evloss=False):
        """ Compute the 6 different phase space factors on the QRPA grid and store the result.

        Args:
            Evloss (bool): Compute phase space for neutrino energy loss rather than beta decay.
                           Only affects result for electron capture.
            contour (famContour): Compute the phase space approximation using a different
                                  contour than the one on which the strength lies. Only affects
                                  result when analyticPsi is called in phase space (i.e. closed contours).
        """

        Zd       = self.sf_metadict[u"Z"]
        A        = self.sf_metadict[u"A"]
        eqrpamax = self.sf_metadict[u"EQRPAmax"]
        hfb_qval = self.sf_metadict[u"HFB_Qval"]
        ft_active= self.sf_metadict[u"ft_active"]
        temper   = self.sf_metadict[u"temper"]

        # Pretend open contours are on real axis, else use complex eqrpa grid
        if not self.contour.closed:
            eqrpa = np.real(self.contour.ctr_z)
        else:
            eqrpa = self.contour.ctr_z
        ps_interval = np.array([self.ps_contour.energy_min, self.ps_contour.energy_max])

        # Compute the phase space integrals at each EQRPA (or approximation thereof)
        ps = pd.DataFrame()
        ps_test = {}
        if self.beta == u'c':
            fac = -1.0
        else:
            fac = +1.0
        w0 = (eqrpamax - eqrpa)/MEC2 + fac # an array
        w0_interval = (eqrpamax - ps_interval)/MEC2 + fac

        # If no first forbidden, don't waste time computing FF phase space
        nps = range(1,7)
        #if all(self.ffk_err): nps = [2] #11/29/23: ignore this condition.

        for i in nps:
            # Q value can be negative for EC at FT
            if hfb_qval < 0 and not (self.beta == u'c' and ft_active):
                ps[u'f'+str(i)] = np.zeros(len(eqrpa))
            else:
                if self.beta == u'c' and ft_active:
                    ps_fct = self.psecFct(i, Zd, A, eqrpamax=eqrpamax, interval=ps_interval, temper=temper,
                            approx=self.contour.closed, Evloss=Evloss)
                elif self.beta == u'c':
                    ps_fct = self.psecFct(i, Zd, A, Bx_shift=False, Evloss=Evloss)
                else:
                    ps_fct = self.psiFct(i, Zd, A, eqrpamax, ps_interval[0], approx=self.contour.closed,
                            Evloss=Evloss)
                psvals = ps_fct(w0)
                ps[u'f'+str(i)] = psvals

                # For debugging
                if self.contour.closed:
                    def mid(x): return (max(np.real(x)) - min(np.real(x)))/2.0 + min(np.real(x))
                    mid_w0 = mid(w0_interval)
                    mid_e  = mid(ps_interval)
                    ps_test[u'f'+str(i)] = {"z_mid":mid_e, "fval":ps_fct(mid_w0)}

        # Store the resulting df
        self.ps_df = ps
        self._ps_test = ps_test

        # Store all the energy data in 1 place
        energy_df = pd.DataFrame()
        energy_df[u'Re(EQRPA)'] = np.real(eqrpa)
        energy_df[u'Im(EQRPA)'] = np.imag(eqrpa)
        energy_df[u'Re(W0)'] = np.real(w0)
        energy_df[u'Im(W0)'] = np.imag(w0)
        self.ps_energy_df = energy_df

        # Store the phase space settings as meta data
        exclude = [u'GA', u'GV']
        exclude_apx = [u'log(pYe)', u'screening', u'psi_glpts']
        if not (self.beta == u'c' and ft_active): exclude += [u'log(pYe)']

        ps_apx_setts = [s for s in self._settings if s not in exclude+exclude_apx]
        ps_setts = {s:self._settings[s] for s in self._settings if s not in exclude}

        # No approx for open contours and terrestrial EC
        if (self.beta == u'c' and not ft_active) or not self.contour.closed:
            ps_setts = {k:None if k in ps_apx_setts else v for k,v in ps_setts.items()}

        self.sf_metadict.update(ps_setts)

    #---------------------------------------------------------------------------
    def prepStrengths(self, zero_neg):
        """ Sorts the raw strengths and crossterms into a dictionary by K
        and multiplies in the necessary dimensionful prefactors.
        Returns (dict): b = { op : {k:cstr_df} }
        """
        def assert_equal(str_a, str_b, label_a):
            """ check crossterms are equal """
            tol_fail = 1e-3; tol_warn = 5e-4
            label_b = u'x'.join(reversed(label_a.split(u'x')))
            diff_cmp = (str_a[label_a].values - str_b[label_b].values)
            diff = {u'Real': abs(np.real(diff_cmp)),
                    u'Imag': abs(np.imag(diff_cmp))}
            for d in diff:
                if any(diff[d] > tol_warn):
                    if any(diff[d] > tol_fail):
                        raise RuntimeError(u"{:} cross-terms {:} do not agree within tolerance.".format(d, label_a))
                    else:
                        print(u"WARNING: {:} cross-terms {:} do not agree very well, but within tolerance.".format(d, label_a))

        dim = len(self.contour.ctr_z)
        # Prefactors (by operator and k)
        Tk  = {u'K=0':1.0, u'K>0':np.sqrt(2)}
        pre = {u'F_K0'   :  Tk[u'K=0'],
               u'GT_K0'  :  Tk[u'K=0'],
               u'GT_K1'  :  Tk[u'K>0'],
               u'P_K0'   :  (-Tk[u'K=0']*HBAR_MEC*MEC2/MN),
               u'R_K0'   : (Tk[u'K=0']*np.sqrt(3.0)/HBAR_MEC),
               u'PS0_K0' :(-1.0*HBAR_MEC*MEC2/MN),
               u'RS0_K0' :  1.0/HBAR_MEC ,
               u'RS1_K0' :(-Tk[u'K=0']/HBAR_MEC),
               u'RS2_K0' : (Tk[u'K=0']/HBAR_MEC),
               u'P_K1'   : (-Tk[u'K>0']*HBAR_MEC*MEC2/MN),
               u'R_K1'   : (Tk[u'K>0']*np.sqrt(3.0)/HBAR_MEC),
               u'RS1_K1' :(-Tk[u'K>0']/HBAR_MEC),
               u'RS2_K1' : (Tk[u'K>0']/HBAR_MEC),
               u'RS2_K2' : (Tk[u'K>0']/HBAR_MEC)}

        #--------------------------------------------------------------
        # Sort raw strength/xterms by operator and k, and add in prefactors.
        #  - Stored as dict:  b = { op : {k:cstr_df} }
        #  - Every operator has k=0,1,2 with 'non-existant' combinations containing
        #    zero strength. This allows us to loop over k=0,1,2 when we compute
        #    the final shape factors C (e.g. k=1,2 for F should be zero)
        #--------------------------------------------------------------
        b = {o : {k : np.zeros(dim) for k in range(3)} for o in \
                shapeFactor.al_ops+shapeFactor.ff_ops+shapeFactor.xts}
        for s in list(self.strengths.values()):
            # Extract useful data for the strength
            genop  = s.genopname # [OP_K#]
            bareop = s.bareop    # [OP]
            k      = s.k
            cstr   = s.cstr_df
            if zero_neg:
                if self.contour.closed:
                    print(u"WARNING: Do not set negative strength to zero for closed contours. Skipping.")
                else:
                    cstr = s.cstr_df_zeroed
            #if self.ffk_err[k] and genop not in shapeFactor.al_ops_raw:
            #    continue
            #11/26 implement new condition to only skip if crossterm components aren't computed. alternately, disable this check altogether
            if genop not in shapeFactor.al_ops_raw and self.ffjk_err(bareop, k):
                print("Skipped because missing crossterms")
                continue
            # Strength
            b[bareop][k] = pre[genop]**2*cstr[u'Strength']

            # Crossterms
            # RS0 x PS0 (k=0)
            if bareop == u'RS0':
                xterm = u'RS0xPS0'
                assert_equal(s.cstr_df, self.strengths[u'PS0_K'+str(k)].cstr_df, xterm)
                b[u'RS0_PS0'][k] = pre[u'RS0_K'+str(k)]*pre[u'PS0_K'+str(k)]*cstr[xterm]
            # R x RS1 (k=0,1)
            if bareop == 'R':
                xterm = 'RxRS1'
                assert_equal(s.cstr_df, self.strengths[u'RS1_K'+str(k)].cstr_df, xterm)
                b[u'R_RS1'][k] = pre[u'R_K'+str(k)]*pre[u'RS1_K'+str(k)]*cstr[xterm]
            # P x RS1 (k=0,1)
            if bareop == u'P':
                xterm = u'PxRS1'
                assert_equal(s.cstr_df, self.strengths[u'RS1_K'+str(k)].cstr_df, xterm)
                b[u'P_RS1'][k] = pre[u'P_K'+str(k)]*pre[u'RS1_K'+str(k)]*cstr[xterm]
            # P x R (k=0,1)
            if bareop == u'R':
                xterm = u'RxP'
                assert_equal(s.cstr_df, self.strengths[u'P_K'+str(k)].cstr_df, xterm)
                b[u'R_P'][k] = pre[u'R_K'+str(k)]*pre[u'P_K'+str(k)]*cstr[xterm]
        #np.savetxt("bb.txt", b[u'PS0'][0]) #TESTSAVE
        return b

    #---------------------------------------------------------------------------
    def calcSfByJ(self, b):
        """ Computes the allowed and first forbidden shape factor from strengths.
        """

        # Terrestrial EC formulas have several sign changes throughout.
        ec      = self.sf_metadict[u"ec"]
        ft_active= self.sf_metadict[u"ft_active"]
        ps      = self.ps_df
        dim     = len(self.contour.ctr_z)
        C_keys  = shapeFactor.C_keys
        all_keys= shapeFactor.al_ops + C_keys
        LAM     = self.sf_metadict[u"LAM"]
        w0max   = self.sf_metadict[u"W0_max"] # Max beta particle energy
        Zd      = self.sf_metadict[u"Z"]      # Z used in SF (B-=Zf, B+=-Zf, EC=Zi)
        rad     = self.sf_metadict[u"Radius"] # Nuclear radius
        mu1     = mu_ke(1)      # Coulomb function (=1)
        gamma1  = gam_ke(1, Zd) # Coulomb function (=sqrt(1 - (alpha*Zd)^2))
        xp      = w0max/3.0 + ec*0.5*ALPHA*Zd/rad # I(1,1,1,1;r) ~ 3/2
        xm      = w0max/3.0 - ec*0.5*ALPHA*Zd/rad # I(1,1,1,1;r) ~ 3/2
        C = {term : {k : np.zeros(dim) for k in range(3)} for term in all_keys}

        for k in range(3):
            #--------------------------------------------------------------
            # Allowed terms (recall b[op][K-DNE] = 0)
            #--------------------------------------------------------------
            C[u'F'][k]  = b[u'F'][k]*ps[u'f2']
            C[u'GT'][k] = LAM**2*b[u'GT'][k]*ps[u'f2']

            #if self.ffk_err[k]: continue #11/29/23: try removing these conditions. may lead to unphysical rates
            #--------------------------------------------------------------
            # J = 0 terms
            #--------------------------------------------------------------
            C[u'J0_R1'][k] = \
                    (-ec*2.0/3.0*LAM**2.0*(xp*b[u'RS0'][k] + b[u'RS0_PS0'][k]))*ps[u'f1']
            #C[u'J0_R1'][k] = \
            #        (-ec*2.0/3.0*LAM**2.0*(xp*b[u'RS0'][k]))*ps[u'f1']
            #C[u'J0_R1'][k] = \
            #        (-ec*2.0/3.0*LAM**2.0*(b[u'RS0_PS0'][k]))*ps[u'f1']
            C[u'J0_R2'][k] = \
                    (LAM**2.0*((xp**2.0 + 1.0/9.0)*b[u'RS0'][k] +\
                    b[u'PS0'][k] + 2.0*xp*b[u'RS0_PS0'][k]))*ps[u'f2']
            #C[u'J0_R2'][k] = LAM**2.0 * (b[u'PS0'][k])*ps[u'f2'] #test: only PS0 contribution.
            #C[u'J0_R2'][k] = (LAM**2.0*(xp**2.0 + 1.0/9.0)*b[u'RS0'][k])*ps[u'f2'] #test: only RS0 
            #C[u'J0_R2'][k] = (LAM**2.0*((xp**2.0 + 1.0/9.0)*b[u'RS0'][k] + b[u'PS0'][k]))*ps[u'f2'] #test: PS0 + RS0 (no crossterms)
            #C[u'J0_R2'][k] = LAM**2.0 * (2.0*xp*b[u'RS0_PS0'][k]) * ps[u'f2']
            #--------------------------------------------------------------
            # J = 1 terms
            #--------------------------------------------------------------
            C[u'J1_R1'][k] = \
                    (-2.0/9.0*(ec*xp*b[u'R'][k] - ec*2.0*LAM**2.0*xm*b[u'RS1'][k]\
                    + LAM*np.sqrt(2.0)*(xp-xm)*b[u'R_RS1'][k]\
                    - ec*np.sqrt(3.0)*b[u'R_P'][k] - LAM*np.sqrt(6.0)*b[u'P_RS1'][k]))*ps[u'f1']
 
            C[u'J1_R2'][k] = \
                    (b[u'P'][k] +\
                    1.0/3.0*xp**2.0*b[u'R'][k] +\
                    2.0/3.0*LAM**2.0*xm**2.0*b[u'RS1'][k] +\
                    1.0/27.0*(b[u'R'][k] + 2.0*LAM**2.0*b[u'RS1'][k] +\
                            ec*2.0*np.sqrt(2.0)*LAM*b[u'R_RS1'][k]) +\
                    np.sqrt(2.0/3.0)*(ec*2.0*LAM*xm*b[u'P_RS1'][k] - np.sqrt(2.0)*xp*b[u'R_P'][k] +\
                            (-ec*2.0/np.sqrt(3.0)*LAM*xm*xp*b[u'R_RS1'][k])))*ps[u'f2']
            if self.beta == u'c' and not ft_active:
                C[u'J1_R2'][k] += -8.0/27.0*(LAM**2.0*b[u'RS1'][k] +\
                        (ec*LAM/np.sqrt(2.0)*b[u'R_RS1'][k]))*ps[u'f1']*ps[u'f3']
            # Regular free coulomb fcts for b+, b-, ecft
            else:
                C[u'J1_R2'][k] += -8.0/27.0*(LAM**2.0*b[u'RS1'][k] +\
                        (ec*LAM/np.sqrt(2.0)*b[u'R_RS1'][k]))*mu1*gamma1*ps[u'f2']
            #C[u'J1_R2'][k] = (np.sqrt(2.0/3.0) * -np.sqrt(2.0)*xp*b[u'R_P'][k] )* ps[u'f2']
            #C[u'J1_R2'][k] = (2.0/3.0*LAM**2.0*xm**2.0*b[u'RS1'][k] + 1.0/27.0*2.0*LAM**2.0*b[u'RS1'][k] -\
            #                   8.0/27.0*(LAM**2.0*b[u'RS1'][k]*mu1*gamma1))*ps[u'f2']
            #C[u'J1_R2'][k] = (1.0/27.0*ec*2.0*np.sqrt(2.0)*LAM*b[u'R_RS1'][k] + np.sqrt(2.0/3.0) *\
            #                   -ec*2.0/np.sqrt(3.0)*LAM*xm*xp*b[u'R_RS1'][k] - 8.0/27.0*ec*LAM/np.sqrt(2.0)*b[u'R_RS1'][k]*mu1*gamma1) * ps[u'f2']
            C[u'J1_R3'][k] = \
                    (4.0/3.0*(np.sqrt(2.0)/3.0*LAM*xp*b[u'R_RS1'][k] +\
                    -ec*2.0/3.0*LAM**2.0*xm*b[u'RS1'][k] +\
                    -np.sqrt(2.0/3.0)*LAM*b[u'P_RS1'][k]))*ps[u'f3']

            C[u'J1_R4'][k] = \
                    (8.0/27.0*LAM**2.0*b[u'RS1'][k])*ps[u'f4']
            
            C[u'J1_R5'][k] = \
                    (1.0/27.0*(2.0*b[u'R'][k] + LAM**2.0*b[u'RS1'][k] +\
                    ec*2.0*np.sqrt(2.0)*LAM*b[u'R_RS1'][k]))*ps[u'f5']
            
            C[u'J1_R6'][k] = \
                    (1.0/27.0*(2.0*b[u'R'][k] + LAM**2.0*b[u'RS1'][k] +\
                    -ec*2.0*np.sqrt(2.0)*LAM*b[u'R_RS1'][k]))*ps[u'f6']
            #--------------------------------------------------------------
            # J = 2 terms
            #--------------------------------------------------------------
            C[u'J2_R5'][k] =  (1.0/9.0*LAM**2.0*b[u'RS2'][k])*ps[u'f5']
            C[u'J2_R6'][k] =  (1.0/9.0*LAM**2.0*b[u'RS2'][k])*ps[u'f6']
        return C

    #---------------------------------------------------------------------------
    def calcSfContributions(self, C):
        """ Calculate contributions by summing over K and operators
        """
        betaout_keys = shapeFactor.betaout_keys
        al_ops       = shapeFactor.al_ops
        CJ0_keys     = shapeFactor.CJ0_keys
        CJ1_keys     = shapeFactor.CJ1_keys
        CJ2_keys     = shapeFactor.CJ2_keys
        C_keys       = shapeFactor.C_keys
        all_keys = al_ops + C_keys

        final_keys = {
            betaout_keys[0]  : {u'keys':all_keys,u'k':[0,1,2] , u'i':0},
            betaout_keys[1]  : {u'keys':al_ops,  u'k':[0,1]   , u'i':1},
            betaout_keys[2]  : {u'keys':[u'GT'], u'k':[0,1]   , u'i':2},
            betaout_keys[3]  : {u'keys':C_keys,  u'k':[0,1,2] , u'i':3},

            betaout_keys[4]  : {u'keys':[u'F'],  u'k':[0]     , u'i':4},
            betaout_keys[5]  : {u'keys':[u'GT'], u'k':[0]     , u'i':5},
            betaout_keys[6]  : {u'keys':[u'GT'], u'k':[1]     , u'i':6},

            betaout_keys[7]  : {u'keys':C_keys,  u'k':[0]     , u'i':7},
            betaout_keys[8]  : {u'keys':C_keys,  u'k':[1]     , u'i':8},
            betaout_keys[9]  : {u'keys':C_keys,  u'k':[2]     , u'i':9},

            betaout_keys[10] : {u'keys':CJ0_keys,u'k':[0,1,2] , u'i':10},
            betaout_keys[11] : {u'keys':CJ1_keys,u'k':[0,1,2] , u'i':11},
            betaout_keys[12] : {u'keys':CJ2_keys,u'k':[0,1,2] , u'i':12},

            betaout_keys[13] : {u'keys':CJ0_keys,u'k':[0]     , u'i':13},
            betaout_keys[14] : {u'keys':CJ1_keys,u'k':[0]     , u'i':14},
            betaout_keys[15] : {u'keys':CJ1_keys,u'k':[1]     , u'i':15},
            betaout_keys[16] : {u'keys':CJ2_keys,u'k':[0]     , u'i':16},
            betaout_keys[17] : {u'keys':CJ2_keys,u'k':[1]     , u'i':17},
            betaout_keys[18] : {u'keys':CJ2_keys,u'k':[2]     , u'i':18}
            }

        # Initialize final df. Stupid way to get a list of the columns in the correct order
        cols_tmp = [(key,d[u'i']) for key,d in list(final_keys.items())]
        cols     = [k[0] for k in sorted(cols_tmp,key= lambda x:x[1])]
        sf_df = pd.DataFrame(columns=cols)

        # Populate the final dataframe by summing over appropriate strengths
        dim = len(self.contour.ctr_z)
        for fk in final_keys:
            tmp = np.zeros(dim, dtype=np.complex_)
            for key2sum in final_keys[fk][u'keys']:
                for k2sum in final_keys[fk][u'k']:
                    tmp += C[key2sum][k2sum]
            sf_df[fk] = tmp

        # Store the result
        self.sf_df = sf_df

    #--------------------------------------------------------------
    @staticmethod
    def findFirstPeak(x, y, bare=False):
        """
        Find the index and value of the first local maximum in a numpy array.
        To avoid false positives (e.g. identifying the trough between two negative peaks)
        this method checks for the first instance of positive slope, positive concatvity
        section leading to a positive peak.

        Note:
            This does not always work! This method needs fixed or eliminated probably.
            Classic problem case is when the first peak is very small and very close to
            an extrememly large peak. In this case the condition that a point is a peak
            if the points to the left and right are both less than the middle will fail
            to detect an actual peak. There is also still trouble not identifying troughs
            between two negative peaks as a positve peak.

            I've looked at various peak finding python functions but have not been able
            to get a reliable result thus far.

        Args:
            x (real ndarray, series): The x-coordinates of the signal.
            y (real ndarray, series): The y-coordinates of the signal.

        Returns:
            int: The index of the x-coordinate data for the first identified peak.
            int: The index of the x-coordinate data for the base of the first peak, by
            travelling backwards from the first peak until sign or concavity changes.
        """

        # Check format of input
        if not isinstance(x, np.ndarray) or not isinstance(y,np.ndarray):
            raise ValueError(u"x and y must be numpy arrays")
        if np.iscomplex(x).all() or np.iscomplex(y).all():
            raise ValueError(u"Inputs must be real valued")
        if len(x) != len(y):
            raise ValueError(u"x and y must be same length")
        else:
            if len(x) < 3:
                print(u"ERROR: Cannot find first peak with less than 3 points.")
                return 0, 0

        # Tolerance for being "negative"
        tol = 1e-4

        # Compute backward finite difference for slope
        dy=np.diff(y,1)
        dx=np.diff(x,1)
        yfirst = (dy/dx)
        # Compute backward finite difference for concavity
        dyfirst = np.diff(yfirst,1)
        dxfirst = np.diff(x[1:],1)
        ysecond = (dyfirst/dxfirst)
        # Use forward difference to assign slope and concavity to initial points
        fd1 = yfirst[0]
        fd2 = (dyfirst[:2]/dx[:2])
        yfirst  = np.insert(yfirst,0,fd1)
        ysecond = np.insert(ysecond,0,fd2)

        # Initialize results
        ind_1stpk = None
        pk_start  = 0
        try:
            # Find positive peaks defined as a point greater than the one to its left AND its right and > 0
            # Exclude 1st point. If bare, just look for 1st peak no matter sign/slop/conc, and return it
            if bare:
                peak_mask = np.r_[False, y[1:] > y[:-1]] & np.r_[y[:-1] > y[1:], True]
                ind_1stpk = np.where(peak_mask)[0][0]
                return ind_1stpk, pk_start
            else:
                pos_peak_mask = np.r_[False, y[1:] > y[:-1]] & np.r_[y[:-1] > y[1:], True] & np.r_[y > tol]
                pos_peak_inds = np.where(pos_peak_mask)[0]

            if len(pos_peak_inds) == 0:
                print(u"WARNING: No positive peaks in strength. Using start of interval.")
                ind_1stpk = 0
            else:
                # Separate peaks at x < 0 and x > 0
                peaks_x_lt0 = [p for p in pos_peak_inds if x[p] < 0]
                peaks_x_ge0 = [p for p in pos_peak_inds if x[p] >= 0]

                if peaks_x_lt0:
                    # Find indices of points with +concavity, +slope, and +value, and also find opposite
                    pos_up_mask = np.r_[y >= tol] & np.r_[yfirst >= 0] & np.r_[ysecond >= 0]
                    pos_up_inds = np.where(pos_up_mask)[0]
                    neg_dn_mask = np.r_[y < tol] | np.r_[yfirst < 0] | np.r_[ysecond < 0]
                    neg_dn_inds = np.where(neg_dn_mask)[0]
                    # Any peaks preceeded by +c+s+v? If so use it. Travel backwards from peak until
                    # the first non +c+s+v point, this defines the start.
                    peak_cands = [i for i in peaks_x_lt0 if i > pos_up_inds[0]]
                    if peak_cands:
                        ind_1stpk = peak_cands[0]
                        neg_b4_pk = [i for i in neg_dn_inds if i < ind_1stpk]
                        if len(neg_b4_pk) > 0: pk_start  = neg_b4_pk[-1] + 1

                if ind_1stpk is None:
                    # If no peaks < 0 satisfying our requirements, use 1st peak > 0.
                    ind_1stpk = peaks_x_ge0[0]
                    # Find indices of negative points
                    neg_mask = np.r_[y < tol]
                    neg_inds = np.where(neg_mask)[0]
                    #Travel backwards from peak until first negative point, this defines the start.
                    neg_b4_pk = [i for i in neg_inds if i < ind_1stpk]
                    if len(neg_b4_pk)>0: pk_start  = neg_b4_pk[-1] + 1

        except Exception as e:
            print(u"WARNING: An error occurred in find 1st peak. Using start of interval instead.")
            print(u"         Error Message: ", e)
            ind_1stpk, pk_start = 0, 0

        return ind_1stpk, pk_start


    #==============================================================#
    #                         Quadrature                           #
    #==============================================================#
    def calcBetaRates(self, emin=None, emax=None, quad=None, zero_neg=False, Evloss=False):
        r"""
        Integrate the shape factor and return the corresponding rates and
        half lives.

        Note:

            The pnfam strength contains the factor of  :math:`-\frac{1}{\pi}`:

            .. math::
                S_{pnfam} = -\frac{1}{\pi} S_{QRPA}

            For closed contours, a complex contour integration is performed. This
            results in needing an extra factor of -1/2:

            .. math::
                \sum_n B_n &= \sum_n Res(S,E_n) = \frac{1}{2 \pi i} \int_c S*dz\\
                         &= \frac{1}{2 \pi i}*(0 + Im(\int_c S*dz))\\
                         &= -\frac{1}{2} * Im(\int_c S_{fam}*dz)

            For open contours, a standard integration is performed, pretending
            that the shape factor lies on the real axis (Note dzdt for open
            contours is np.ones(nr_points)):

            .. math::
              \frac{dB}{dw}  &= -\frac{1}{\pi} Im(S)\\
              B_{total} &= \int Im(S_{fam}) dE

        Args:
            emin (float): Custom mininum energy for the integration limits on open contours.
                Ignored for closed contours. (default None).
            emax (float): Custom maximum energy for the integration limits on open contours.
                Ignored for closed contours. (default None).
            quad (str): Quadrature type for the integration. Implemented types: TRAP, GAUSS, SIMPSON.

        Returns:
            dataframe: The rate and half-life contributions to be written to beta.out.
        """
        # Setup/checks
        if quad is None:
            quad = self.contour.quadrature

        if self.sf_df_imag is None:
            raise RuntimeError(u"calcShapeFactor must be run before calcBetaRates.")

        beta_df = self.sf_df_real + self.sf_df_imag*1j
        beta_df.drop([u'Re(EQRPA)', u'Im(EQRPA)'], axis=1, inplace=True)

        headers = [u'Rate(s^-1)', u'Half-Life(s)', u'Contribution']
        betaout = defaultdict(list)

        for h in beta_df:
            if h in shapeFactor.betaout_keys:
                ctr_int = self.complex_quadrature(quad, self.contour, beta_df[h], emin, emax)
                if self.contour.closed:
                    # Need Re[cint S], df is -1/pi S. Im[ i * cint] to avoid -0 default val
                    rate = (np.log(2)/KAPPA*np.imag(-1j*np.pi*ctr_int))
                else:
                    rate = (np.log(2)/KAPPA*np.imag(ctr_int)) # Need -1/pi Im[int S]
                # Suppress warning if divide by zero, but keep the inf result
                with np.errstate(divide=u'ignore'):
                    halflife = (np.log(2)/rate)

                betaout[headers[0]].append(rate)
                betaout[headers[1]].append(halflife)
                betaout[headers[2]].append(h)

        # Make dataframe, order columns, set contributation as index
        beta_df = pd.DataFrame(betaout)
        beta_df = beta_df[headers]
        beta_df.set_index(headers[2], inplace=True)

        self.sf_metadict.update({u'quadratr':quad, u'FAM_ctr':self.contour.name_and_int})
        return beta_df

    #--------------------------------------------------------------
    @staticmethod
    def zeroNegRatesBetaOut(beta_df):
        """ Given a beta.out dataframe, zero negative rate contributions and resum
            This should make no change for open contours, which have the strength
            functions themselves set to zero where negative.
        """
        beta_df = beta_df.copy()
        headers = [u'Rate(s^-1)', u'Half-Life(s)', u'Contribution']
        beta_totals      = shapeFactor.beta_totals
        beta_al_contribs = shapeFactor.beta_al_contribs
        beta_ffK_totals  = shapeFactor.beta_ffK_totals
        beta_ffJ_totals  = shapeFactor.beta_ffJ_totals
        beta_ff_contribs = shapeFactor.beta_ff_contribs
        totals = {beta_totals[0]    : beta_al_contribs+beta_ff_contribs,
                  beta_totals[1]    : beta_al_contribs,
                  beta_totals[2]    : [c for c in beta_al_contribs if u'GT' in c],
                  beta_totals[3]    : beta_ff_contribs,
                  beta_ffK_totals[0]: [c for c in beta_ff_contribs if c[-2]==u'0'],
                  beta_ffK_totals[1]: [c for c in beta_ff_contribs if c[-2]==u'1'],
                  beta_ffK_totals[2]: [c for c in beta_ff_contribs if c[-2]==u'2'],
                  beta_ffJ_totals[0]: [c for c in beta_ff_contribs if c[-4]==u'0'],
                  beta_ffJ_totals[1]: [c for c in beta_ff_contribs if c[-4]==u'1'],
                  beta_ffJ_totals[2]: [c for c in beta_ff_contribs if c[-4]==u'2']}

        # Zero negative rates (allow beta_df to have both Rate and Evloss-Rate)
        rtcols = [h for h in beta_df if headers[0] in h]
        hlcols = [h for h in beta_df if headers[1] in h]
        beta_df[beta_df.loc[:, rtcols] < 0] = 0
        with np.errstate(divide=u'ignore'):
            beta_df.loc[:, hlcols] = np.log(2)/beta_df.loc[:, rtcols].values

        # Resum totals
        for hrt, hhl in zip(rtcols, hlcols):
            for h in totals:
                new_rt = 0
                for c in totals[h]:
                    new_rt += beta_df.loc[c,hrt]
                beta_df.loc[h, hrt] = new_rt
                with np.errstate(divide=u'ignore'):
                    beta_df.loc[h, hhl] = np.log(2)/new_rt

        return beta_df

    #--------------------------------------------------------------
    @staticmethod
    def recalcHalflives(beta_df):
        beta_df = beta_df.copy()
        headers = [u'Rate(s^-1)', u'Half-Life(s)', u'Contribution']
        rtcols = [h for h in beta_df if headers[0] in h]
        hlcols = [h for h in beta_df if headers[1] in h]
        with np.errstate(divide=u'ignore'):
            beta_df.loc[:, hlcols] = np.log(2)/beta_df.loc[:, rtcols].values
        return beta_df

    #--------------------------------------------------------------
    def calcTotalStr(self, emin=None, emax=None, quad=None, zero_neg=False):
        # Setup/checks
        if quad is None:
            quad = self.contour.quadrature

        # Include squared prefactors
        # This includes intrinsic to lab frame, so theta_k (e.g. B[GTK1] = 2 GTK1)
        # This includes dimensional prefactors so operators are dimensionless (PRC 90 024308 (2014))
        # This does NOT include GA/LAM (Ikeda sum rule does not include GA...)
        b = self.prepStrengths(zero_neg)

        # Integrate strengths of operators we calculated
        bout = defaultdict(list)
        headers = [u'B_total', u'Unphysical', u'Operator']
        for s in list(self.strengths.values()):
            cint = self.complex_quadrature(quad, self.contour, b[s.bareop][s.k], emin, emax)
            # Closed contours:
            #     B = RES[S] = Re[(cint S)/(2pi i)]
            #     Include factor of -pi since we integrated (-1/pi S), not S
            #     Include factor of i so B = Im[cint(S)*(-pi*i)] to match open contours
            # Open contours:
            #     On real axis, dB/dw = -1/pi Im(S), so B = Im[int(-1/pi S)]
            if self.contour.closed:
                cint *= (-np.pi*1j)
            btot = np.imag(cint)
            unphys = np.real(cint)

            bout[headers[0]].append(btot)
            bout[headers[1]].append(unphys)
            bout[headers[2]].append(s.opname)

        # Make dataframe, order columns, set contributation as index
        b_df = pd.DataFrame(bout)
        b_df = b_df[headers]
        b_df.set_index(u'Operator', inplace=True)

        return b_df

    #--------------------------------------------------------------
    def integratePhaseSpace(self, emin=None, emax=None, quad=None):
        if quad is None:
            quad = self.contour.quadrature

        headers = [u'Re[Int_f]', u'Im[Int_f]', u'Factor']
        pout = defaultdict(list)

        def ftest(x, x0): return res/(x - x0)

        for fps in self.ps_df:
            try:
                pole = self._ps_test[fps]["z_mid"]
                fval = self._ps_test[fps]["fval"]
                res = 1.0 + 2.0j
                test = ftest(self.contour.ctr_z,pole)
            except KeyError:
                res, test, fval = 0.0, 1.0, 1.0
            if abs(fval) < 1e-20:
                cint = self.complex_quadrature(quad, self.contour, self.ps_df[fps], emin, emax)
            else:
                cint = self.complex_quadrature(quad, self.contour, test*self.ps_df[fps], emin, emax)
                cint = cint/fval - res # Should be zero if no other poles...
            pout[headers[0]].append(np.real(cint))
            pout[headers[1]].append(np.imag(cint))
            pout[headers[2]].append(fps)

        p_df = pd.DataFrame(pout)
        p_df = p_df[headers]
        p_df.set_index(u'Factor', inplace=True)

        self.sf_metadict.update({u'quadratr':quad, u'FAM_ctr':self.contour.name_and_int})
        return p_df

    #--------------------------------------------------------------
    @staticmethod
    def complex_quadrature(quad, contour, y, xmin=None, xmax=None):
        r"""
            Regular integration or complex contour integration of y-data along a
            given contour. CCI is performed if the contour is closed.
        """
        # Setup/checks
        if quad not in [u'TRAP', u'GAUSS', u'SIMPSON']:
            raise ValueError(u"Invalid quadrature requested.")

        # Real or complex quadrature
        if contour.closed:
            fac  = 1.0/(2.0*np.pi*1j) # [.5/pi Y' - i .5/pi Y]    and Y=-Y/pi
            dzdt = contour.ctr_dzdt
            dt   = contour.theta
        else:
            fac  = 1.0
            dzdt = 1.0
            dt   = np.real(contour.ctr_z)
            # Optionally adjust the domain
            if xmin is None: xmin = min(dt)
            if xmax is None: xmax = max(dt)
            mask = np.r_[dt>=xmin] & np.r_[dt<=xmax]
            dt = dt[mask]
            y = y[mask]

        if quad == u'GAUSS':
            if not contour.use_gauleg:
                raise RuntimeError(u"Integration points do not lie on a gauss-legendre grid.")
            glwts= contour.glwts
            cint = fac*sum((y*dzdt)*glwts)
        elif quad == u'TRAP':
            cint = fac*trapz(y*dzdt, dt)
        elif quad == u'SIMPSON':
            cint = fac*simps(y*dzdt, dt)

        return cint

    #--------------------------------------------------------------
    def calcCumulativeStr(self, str_df=None, emin=None, emax=None, rate=True):
        r"""
        Cumulative integration of the strength for open contours using trapezoidal rule.

        Args:
            str_df (dataframe):
            emin (float): Custom mininum energy for the integration limits (default None).
            emax (float): Custom maximum energy for the integration limits (default None).
            rate (bool): Include the dimensionful prefactor :math:`\frac{\ln2}{\kappa}`.
                One might not want to include this if calculating cumulative bare strength,
                but would want to include this if calculating cumulative rate (beta intensity).
                (default True).

        Returns:
            DataFrame
        """

        if self.contour.closed:
            raise RuntimeError(u"Cumulative strength should not be calculated for closed contour data")

        # Get DF in desired range
        if str_df is None: str_df = self.sf_df_imag
        if emin   is None: emin = min(str_df[u'Re(EQRPA)'])
        if emax   is None: emax = max(str_df[u'Re(EQRPA)'])
        i_min = str_df[str_df[u'Re(EQRPA)'] >= emin].index[0]
        i_max = str_df[str_df[u'Re(EQRPA)'] <= emax].index[-1]
        str_df = str_df.loc[i_min:i_max]

        cumstr_out = pd.DataFrame()
        cumstr_out[u'Re(EQRPA)'] = str_df[u'Re(EQRPA)']
        if rate:
            fac = (np.log(2)/KAPPA)
        else:
            fac = 1.0

        for h in str_df:
            if h in shapeFactor.betaout_keys or h == 'Total-SD':
                cumulative = cumtrapz(str_df[h], str_df[u'Re(EQRPA)'], initial=0)
                cumstr_out[h] = cumulative*fac

        self.sf_metadict.update({u'quadratr':u'TRAP', u'FAM_ctr':self.contour.name_and_int})
        return cumstr_out

    #--------------------------------------------------------------
    def calcTotalGT(self, emin=None, emax=None, zero_neg=False):
        """
        Calculate the total GamowTeller strength as GT(K=0) + 2*GT(K=1). Note
        This does NOT include the appropriate factor of lambda^2.

        Args:
            emin (float): Custom mininum energy for the integration limits (default None).
            emax (float): Custom maximum energy for the integration limits (default None).
            zero_neg (bool): Set negative strength to zero before totalling (default False).

        Returns:
            DataFrame
        """

        genop_list = list(self.strengths.keys())
        if u'GT_K0' not in genop_list or u'GT_K1' not in genop_list:
            raise RuntimeError(u"Missing GT strength for calcTotalGT")

        gt_df = pd.DataFrame(columns=[u'Re(EQRPA)', u'Im(EQRPA)'])
        gt_df[u'Re(EQRPA)'] = np.real(self.contour.ctr_z)
        gt_df[u'Im(EQRPA)'] = np.imag(self.contour.ctr_z)

        gt0 = np.imag(self.strengths[u'GT_K0'].cstr_df[u'Strength'].values)
        gt1 = np.imag(self.strengths[u'GT_K1'].cstr_df[u'Strength'].values)

        gt_df[u'Total-GT'] = gt0 + 2.0*gt1
        if zero_neg:
            gt_df[u'Total-GT'].iloc[gt_df.index[gt_df[u'Total-GT'] < 0]] = 0

        gtc = self.calcCumulativeStr(gt_df, emin=emin, emax=emax, rate=False)
        gt_df[u'Cumulative-GT'] = gtc[u'Total-GT']

        return gt_df

    #--------------------------------------------------------------
    def calcTotalSD(self, emin=None, emax=None, zero_neg_str=False):
        """
        Calculate the total spin dipole strength.

        Args:
            emin (float): Custom mininum energy for the integration limits (default None).
            emax (float): Custom maximum energy for the integration limits (default None).
            zero_neg_str (bool): Set negative strength to zero before totalling (default False).

        Returns:
            DataFrame
        """

        genops = self.strengths.keys()
        if all([not {u'RS0_K'+k, u'RS1_K'+k, u'RS2_K'+k} <= genops for k in (u'0', u'1', u'2')]):
            raise RuntimeError(u"Missing SD strength for calcTotalSD")

        sd_df = pd.DataFrame(columns=[u'Re(EQRPA)', u'Im(EQRPA)'])
        sd_df[u'Re(EQRPA)'] = np.real(self.contour.ctr_z)
        sd_df[u'Im(EQRPA)'] = np.imag(self.contour.ctr_z)
        total_SD = np.array((0.0,) * len(sd_df))

        for j in (0, 1, 2):
            for k in range(0, j+1):
                item = u'RS'+str(j)+'_K'+str(k)
                factor = 1.0 if k == 0 else 2.0
                if item in genops:
                    total_SD += factor * np.imag(self.strengths[item].cstr_df[u'Strength'].values)

        if zero_neg_str:
            total_SD.clip(min=0)
        sd_df[u'Total-SD'] = total_SD

        sdc = self.calcCumulativeStr(sd_df, emin=emin, emax=emax, rate=False)
        sd_df[u'Cumulative-SD'] = sdc[u'Total-SD']

        return sd_df

    #==============================================================#
    #                       I/O Methods                            #
    #==============================================================#
    #--------------------------------------------------------------
    def sfSummaryString(self, sd_in):
        """
        Construct a formatted header string with some key info for output files.

        Args:
            sd_in (dict): Dictionary with key info parameter names/values.

        Returns:
            str
        """

        # Proper formatting for type
        sd = {}
        for k in sd_in:
            if sd_in[k] is None:
                sd[k] = (u'{:<10}', u'N/A')
            elif isinstance(sd_in[k], (float, np.floating)):
                sd[k] = (u'{:<10.5f}', sd_in[k])
            elif isinstance(sd_in[k], (int, np.integer)):
                sd[k] = (u'{:<10d}', sd_in[k])
            else:
                # Smallest multiple of 10 greater than object len
                try:
                    slen = max(10,len(sd_in[k]))
                except:
                    slen = 10
                if slen > 10: slen = 26+10
                sd[k] = (u'{:<'+str(slen)+'}', sd_in[k])

        # Adjust some names for phase space variables
        if sd_in[u'psi_approx'] == u'RATINT':
            psi_line = [u'ratint_pts']
        elif sd_in[u'psi_approx'] == u'POLYFIT':
            psi_line = [u'polyfit_pts', u'polyfit_ord']
            sd[u'polyfit_pts'] = sd.pop(u'polynomial_fit_pts')
            sd[u'polyfit_ord'] = sd.pop(u'polynomial_order')
        else:
            psi_line = []
        if u'log(pYe)' in sd:
            psi_line += ['log(pYe)']

        # Define variables in each line
        keys_by_line = [[u'FAM_ctr',    u'temper'],
                        [u'beta_type',  u'quadratr',   u'Half_Width', u'screening'],
                        [u'psi_approx', u'psi_glpts']+psi_line,
                        [u'Zi',         u'A',          u'Zf',         u'HFB_Qval'],
                        [u'FAM_Qval',   u'EQRPAmax',   u'E_1stPeak',  u'|gA|/gV'],
                        [u'gA',         u'gV',         u'M_nucleon',  u'alpha*Z'],
                        [u'Radius',     u'alpha*Z/2R', u'W0_max',     u'W0*R']
                        ]

        # contruct the string
        date = str(datetime.datetime.now().replace(microsecond=0))[:-3]
        header = u'# PynFAM code version: {:}\n'.format(__version__)+\
                 u'# Run Date: {:}\n'.format(date)+\
                 u'# Summary Data:\n'
        body = u''
        for keys in keys_by_line:
            string = u', '.join([u'{:<11} = '+sd[k][0] for k in keys])
            vals   = [val for k in keys for val in (k, sd[k][1])]
            body += (u'#   '+string+u'\n').format(*vals)
        footer = u'#'

        return header+body+footer

    #--------------------------------------------------------------
    def writeOutput(self, df, title, fname, dest=u'./'):
        """
        Write a data frame with summary data at top to file
        (for beta.out, shapefactor.out, phasespace.out, ...).

        Args:
            df (dataframe): Data to be written
            title (str): Title for the contents, written on the first line of the output file.
            fname (str): Filename for the output file.
            dest (str): Destination for the output file (default './')
        """
        pd_string = df.to_string(header=True, index=True, col_space=3,
                                 float_format=lambda x: u'{:25.16e}'.format(x), index_names=False)

        with open(os.path.join(dest, fname),u'w') as f:
            f.write(title+u'\n')
            f.write(self.sf_metastr+u'\n')
            f.write(pd_string+u'\n')
