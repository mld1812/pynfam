# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
# -------------- Utilities -----------------
import numpy as np
import os
# ----------- Relative Imports -------------
from .fortran_utils          import fortProcess
from ..outputs.ls_logger     import hfbLogger
from ..outputs.hfbtho_parser import hfbthoParser
from ..config                import DEFAULTS
from ..config                import DMASS_NH, MEC2
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS hfbthoRun                                   #
#===============================================================================#
class hfbthoRun(fortProcess):
    """
    An hfbthoRun is a fortProcess subclass that contains the methods to
    run the HFBTHO fortran program and handle the output.

    Args:
        paths (pynfamPaths): The paths for the pynfam calculation.
        beta_type (str): Beta decay type (default None).

    Attributes:
        log (hfbLogger): Logger object for logging key outputs.
        core_rundir (str): The run directory label of the even core for blocking.
        beta (str): Beta decay type, for determining Q value.
        _sd (dict): Internal initialization for soln_dict attribute.
        _soln (hfbthoParser): Internal initialization for soln attribute.
    """

    file_exe=u"hfbtho_main"
    """ file_exe (str): File name of executable.
    """
    file_nml=u"hfbtho_NAMELIST.dat"
    """file_nml (str): File name of namelist.
    """
    file_txt=u"thoout.dat"
    """file_txt (str): File name of textfile output.
    """
    file_bin=u"hfbtho_output.hel"
    """file_bin (str): File name of binary output.
    """
    #file_fam=u"solution.hfb"
    #"""file_fam (str): File name of binary output for fam.
    #"""
    file_log=u"logfile_mini.dat"
    """file_log (str): File name of mini logfile.
    """
    file_tar=u"hfb_meta_solns.tar"
    """file_log (str): File name of mini logfile.
    """

    def __init__(self, paths, beta_type=None):
        # Inheret/override attributes of fortran process(includes output):
        exe = hfbthoRun.file_exe
        fortProcess.__init__(self, paths, exe, hfbthoRun.file_nml, DEFAULTS[u'hfb'])

        # Processing
        self._sd = {u'Z-Blk':u'Error', u'N-Blk':u'Error'}
        self.output = None
        self.soln_err = False
        self.soln_err_msg = u''
        self.core_rundir = None

        # Allow specifying beta+ or beta-. If None then both are stored in log.
        self.beta = beta_type

    @property
    def soln_dict(self):
        """ soln_dict (dict): Contains key results parsed from outputs, if available.
        Always contains label and nucleus.
        """
        sd = self._sd
        sd[u'Label'] = self.label
        sd[u'Z']     = self.nucleus[1]
        #sd[u'Z-Blk'] = self.blocking_str[1]
        sd[u'N']     = self.nucleus[0]
        #sd[u'N-Blk'] = self.blocking_str[0]
        sd[u'b2']    = self.basisDef
        if self.ft_active:
            sd[u'Temp'] = self.temperature
        return sd

    @property
    def nucleus(self):
        """ nucleus (tuple): Neutrons, protons, mass (N, Z, A).
        """
        Z = self.nml[u'HFBTHO_GENERAL'][u'proton_number'] + \
            np.sign(self.nml['HFBTHO_BLOCKING'][u'proton_blocking'][0])
        N = self.nml[u'HFBTHO_GENERAL'][u'neutron_number'] + \
            np.sign(self.nml[u'HFBTHO_BLOCKING'][u'neutron_blocking'][0])
        A = N+Z
        return (N, Z, A)

    @property
    def blocking(self):
        """ blocking (tuple): (n-blocking list, p-blocking list).
        """
        return (self.nml[u'HFBTHO_BLOCKING'][u'neutron_blocking'],
                self.nml[u'HFBTHO_BLOCKING'][u'proton_blocking'])

    @property
    def constraints(self):
        """ constraints (list): Expectation values.
        """
        return self.nml[u'HFBTHO_CONSTRAINTS'][u'expectation_values']

    @property
    def activeConstraints(self):
        """ activeConstraints (list): Lambda active.
        """
        return self.nml[u'HFBTHO_CONSTRAINTS'][u'lambda_active']

    @property
    def basisDef(self):
        """ basisDef (float): Basis deformation.
        """
        return self.nml[u'HFBTHO_GENERAL'][u'basis_deformation']

    @property
    def even(self):
        """ even (bool): True if nucleus is even-even.
        """
        return (self.nucleus[0]%2==0 and self.nucleus[1]%2==0)

    @property
    def blocking_str(self):
        """ blocking_str (tuple): Contains a string with the asymptotic quantum
        numbers of the blocked state {2*omega}{parity}[nn=nz+2*nr+nl, nz, nl].
        We keep the sign on omega to indicate the direction of blocking.
        """
        return self.get_blocking_str(self.blocking)

    @property
    def ft_active(self):
        """ ft_active (bool): True for finite-temperature calculations.
            NB! HFBTHOv2 paper is wrong. Temperature=0 is always used if
            set_temperature=F. If set_temperature=T and temperature<1e-10,
            set_temperature=F is imposed.
        """
        ftbool = self.nml[u'HFBTHO_TEMPERATURE'][u'set_temperature']
        temper = self.nml[u'HFBTHO_TEMPERATURE'][u'temperature']
        return (ftbool and temper > 1e-10)

    @property
    def temperature(self):
        """ ft_active (bool): True for finite-temperature calculations.
            NB! HFBTHOv2 paper is wrong. Temperature=0 is always used if
            set_temperature=F. If set_temperature=T and temperature<1e-10,
            set_temperature=F is imposed.
        """
        if self.ft_active:
            return self.nml[u'HFBTHO_TEMPERATURE'][u'temperature']
        else:
            return 0.0

    @staticmethod
    def get_blocking_str(blocking_in):
        """ blocking_str (tuple): Contains a string with the asymptotic quantum
        numbers of the blocked state {2*omega}{parity}[nn=nz+2*nr+nl, nz, nl].
        We keep the sign on omega to indicate the direction of blocking.
        """
        bs = [u"--",u"--"]
        try:
            for i, blocking in enumerate(blocking_in):
                if blocking[0]:
                    omega = u"{:>+d}".format(blocking[0]) # Keep sign, indicates direction to block
                    par = u"?" # In case we are doing a scan (which we should never do) this is not defined.
                    if np.sign(blocking[1]) > 0:
                        par = u"+"
                    else:
                        par = u"-"
                    qns = u"[{:>d},{:>d},{:>d}]".format(*blocking[2:])
                    bs[i] = omega+par+qns
        except:
            bs = 2*[u"Error"]
        return tuple(bs)

    #-----------------------------------------------------------------------
    def updateOutput(self, output, get_soln=True):
        """
        Wrapper to bundle setting output and parsing it for solution values.

        If get_soln is True, getSoln(output) is called and the output is
        stored in soln_dict, soln_err, soln_err_msg attributes.

        Args:
            output (list of str, or str): A list of the hfbtho output contents,
                str is assumed to be the path to the output file.
            get_soln (bool): Parse the output and populate the solution
                attributes (default True).
        """
        self.output = output
        if get_soln:
            sd, err, emsg = self.getSoln(output)
            self.soln_dict.update(sd)
            self.soln_err = err
            self.soln_err_msg = emsg

    #-----------------------------------------------------------------------
    def setTemperature(self, temp):
        """
        Update the namelist values for a given nucleus.

        Args:
            temp (float): Temperature in MeV
        """
        if not temp or temp < 1e-10:
            set_temp = False
            temp = 0.0
        else:
            set_temp = True

        nuc_inputs = {u'HFBTHO_TEMPERATURE':
                        {u'set_temperature': set_temp,
                         u'temperature':temp}
                     }
        self.adjustNml(nuc_inputs)

    #-----------------------------------------------------------------------
    def setNucleus(self, N, Z):
        """
        Update the namelist values for a given nucleus.

        Args:
            N (int): Number neutrons.
            Z (int): Number protons.
        """

        nuc_inputs = {u'HFBTHO_GENERAL':
                        {u'proton_number':Z,u'neutron_number':N}
                     }
        self.adjustNml(nuc_inputs)

    #-----------------------------------------------------------------------
    def setConstraint(self, lambda_val, expect_val, kickoff=1):
        """
        Update the namelist values for a given constraint.

        Args:
            lambda_val (int): Which constraint to set.
            expect_val (float): Value to set.
            kickoff (int): Constraint type (default 1).
        """

        active, expect_vals = [None]*8, [None]*8
        # Active is auto set to 1, can be set to -1 or 0 with override
        active[lambda_val] = kickoff
        expect_vals[lambda_val] = expect_val
        constraint_inputs = {u'HFBTHO_CONSTRAINTS':
                                {u'lambda_active': active,
                                 u'expectation_values':expect_vals}
                            }
        self.adjustNml(constraint_inputs)

    #-----------------------------------------------------------------------
    def setBeta2s(self, B2):
        """
        Update the namelist values for basis deformation and initial woods-
        saxxon deformation with the same beta2 value.

        Args:
            B2 (float): Deformation beta2.
        """

        basis_inputs = {u'HFBTHO_GENERAL':
                           {u'basis_deformation':B2 },
                        u'HFBTHO_INITIAL':
                           {u'beta2_deformation':B2 }
                        }
        self.adjustNml(basis_inputs)

    #-----------------------------------------------------------------------
    def calculatedQ2(self):
        """
        First order mass dependent calculation of Q2 from beta2.

        Returns:
            float
        """

        A  = self.nucleus[2]
        B2 = self.basisDef
        Q2 = float(B2)*((float(A))**(5.0/3.0))*np.sqrt(5.0/np.pi)*0.01
        return Q2

    #-----------------------------------------------------------------------
    def setRestartFile(self, restart):
        """
        Update the namelist restart from file option (boolean input accepted).

        Args:
            restart (int, bool): Bool is converted to int for HFBTHO.
        """

        if isinstance(restart, bool):
            if restart:
                restart = -1
            else:
                restart = 1
        iter_inputs = {u'HFBTHO_ITERATIONS':
                           {u'restart_file':restart },
                        }
        self.adjustNml(iter_inputs)

    #-----------------------------------------------------------------------
    def setBlocking(self, n, p, key=None):
        """
        Update the namelist values with blocking info.

        * key=None: n,p are lists of 5 quantum numbers for candidate.
        * key=0: n,p are the first quantum number :math:`(2*\Omega)`
            (!=0 turns blocking on, sign determines :math:`A\pm1`).
        * key=1: n,p are the last 4 quantum numbers :math:`(\Pi, N, n_z, n_r)`.

        Args:
            n (list of int, int): See above notes.
            p (list of int, int): See above notes.
            key (int): (default None): See above notes.
        """

        err = None
        if key is None:    # Set all 5 blocking params
            if len(p)!=5 or len(n)!= 5: err = 1
            p_list = p
            n_list = n
        elif key == 0: # Set only 1st blocking (2*Omega, on/off, A+/-1)
            if not isinstance(p,int) or not isinstance(n,int): err = 1
            p_list = [p]+[None]*4
            n_list = [n]+[None]*4
        elif key == 1: # Set only last 4 blocking (parity, N, n_z, n_r)
            if len(p)!=4 or len(n)!=4: err = 1
            p_list = [None]+p
            n_list = [None]+n
        else:
            raise ValueError(u"Invalid key for method setBlocking.")
        if err:
            raise ValueError(u"Inputs must be appropriate length/type for provided key in method setBlocking.")

        blocking_inputs = {u'HFBTHO_BLOCKING':{u'proton_blocking':p_list,
                                              u'neutron_blocking':n_list}
                          }
        self.adjustNml(blocking_inputs)

    #-----------------------------------------------------------------------
    def increment_nucleus(self, dripline, n=1, direction=u'N'):
        """
        Increase proton or neutron number according to dripline input parameter.
        Note dripline=0 is valid, in which case this method does not change anything.

        Args:
            dripline (int): Dripline mode value, as in the input file.
            n (int): Integer multiple, increments n*abs(dripline) (default 1).
            direction (str): Increase neutron ('N') or proton ('P') number (default 'N').
        """

        if not dripline: return

        if direction == u'N':
            N_prev = self.nucleus[0]
            hfbtho_number = u'neutron_number'
            hfbtho_blocking = u'neutron_blocking'
        else:
            N_prev = self.nucleus[1]
            hfbtho_number = u'proton_number'
            hfbtho_blocking = u'proton_blocking'

        block = np.sign(dripline) # (+1, -1, 0)
        increment = n*abs(dripline) # n*(0, 1, 2)
        N_next = N_prev + increment

       # HFBTHO nucleus must always be even, hfb_obj.nucleus includes +/-1 for
        # blocking and therefore could be odd, so adjust accordingly.
        if N_next%2 == 0:
            N_hfb = N_next
            b_hfb = 0
        else:
            N_hfb = N_next - block
            b_hfb = block

        new_params = {
            hfbtho_number: N_hfb,
            hfbtho_blocking: [b_hfb,0,0,0,0]
            }
        self.setNmlParam(new_params)

    #-----------------------------------------------------------------------
    def getBlockingResult(self, output=None):
        """ Parse the asymptotic quantum numbers of the blocked state AFTER
        self consistent loop convergence. This may be different than the quantum
        numbers requested in the namelist.
        """

        no_blo = [0,0,0,0,0]
        blocking = [list(no_blo), list(no_blo)]
        #if not self.blocking[0][0] and not self.blocking[1][0]:
        #    return blocking

        if output is None: output = self.output
        soln = hfbthoParser(output)
        pn_data = soln.getQPData()
        if pn_data is None:
            return None

        blo_data = [pn_data[u'n'][u'blo_qp'], pn_data[u'p'][u'blo_qp']]
        # {2*omega}{parity}[nn=nz+2*nr+nl, nz, nl].
        for i, bd in enumerate(blo_data):
            if bd is not None:
                blocking[i][0] = np.sign(self.blocking[i][0])*abs(bd[u'omega'])
                blocking[i][1] = bd[u'par']
                blocking[i][2] = bd[u'N']
                blocking[i][3] = bd[u'nz']
                blocking[i][4] = bd[u'nl']

        return tuple(blocking)

    #-----------------------------------------------------------------------
    def getSpinPar(self, output=None):
        r"""
        Construct a string for the nucleus' spin and parity.

        Uses the Gallagher-Moskowski Rule for determining spin-parity in
        odd-odd nuclei (Gallagher, Moszkowski, Phys. Rev., 111, 5, 1958).
        The result depends on the spin quantum numbers of the 2 blocked states:

            * spins aligned : :math:`K = K_1 + K_2`
            * spins anti-aligned: :math:`K = K_1 + \overline{K_2} = K_1 - K_2`

        Args:
            output (list of str, str): Input for hfbthoParser.

        Returns:
            str: 'KP', K is in {.1f} form for half interger, P is '+' or '-'.
        """

        if output is None: output = self.output

        soln = hfbthoParser(output)

        pn_data = soln.getQPData()
        if pn_data is None:
            return "Error"

        # Even-Even are 0+
        Ki = 0; Pi = 1; Pi_str = u'+'

        pblo = pn_data[u'p'][u'blo_qp']
        nblo = pn_data[u'n'][u'blo_qp']

        # Odd-Odd
        if pblo is not None and nblo is not None:
            if pblo[u'ns'] == nblo[u'ns']:                # Spins aligned
                Ki = pblo[u'omega'] + nblo[u'omega']      # K1 + K2
            else:                                         # Spins anti-aligned
                Ki = abs(pblo[u'omega'] - nblo[u'omega']) # K1 + \bar{K2}
            Pi = pblo[u'par']*nblo[u'par']
        # Proton-Odd
        elif pblo is not None:
            Ki = pblo[u'omega']
            Pi = pblo[u'par']
        # Neutron-Odd
        elif nblo is not None:
            Ki = nblo[u'omega']
            Pi = nblo[u'par']

        # Construct the string for the logfile
        if Pi < 0: Pi_str = u'-'

        return u'{:.1f}{:}'.format(Ki/2.0, Pi_str)

    #-----------------------------------------------------------------------
    def getEqrpaMax(self, beta, output=None):
        r"""
        Calculate the maximum QRPA energy relevant for spontaneous beta decay.

        Value is determined from the Q-value approximation in Ref.
        (Engel et al., Phys. Rev. C 60, 014302 (1999)).

        .. math::
            Q_{\beta^-} = \Delta M_{n-H} + B_{g.s.}(Z,N) - B_{g.s.}(Z+1,N-1) \\
            B_{g.s.}(Z+1,N-1) \approx B_{g.s.}(Z,N) + \lambda_p - \lambda_n + E_{2qp,lowest} \\
            Q_{\beta^-} \approx [\Delta M_{n-H} + \lambda_n - \lambda_p] - E_{2qp,lowest}

        Note:
            The max QRPA energy changes in a magic nucleus (i.e. when pairing
            collapses). In this case we should use the correct s.p. energies rather than
            the arbtrary q.p. energies (recall when pairing collapses, N is a good
            quantum number and the HFBTHO constraint 'ala' becomes arbitrary. We should
            really use 'alast', which chooses the s.p. fermi energy if pairing collapse
            is detected). However:

               1. A change in the lambdas will cancel in Q.
               2. A change in the lambdas will cause the EQRPA_max and E_gs change in
                  the same way, so the beta decay energy interval stays the same.
               3. The beta-decay phase space integrals (at least for beta-minus...haven't
                  checked beta-plus) are unchanged with a change of lambdas.
               4. The QRPA phonon spectrum shifts in the same way as E_gs and
                  EQRPA_max, so the relative location of peaks does not change.

            As a result, there's no need to make a correction, and we can continue to use
            the QP fermi energies. The strength function will look exactly the same, just
            shifted, but the EQRPA_max and E_gs will also be shifted in the same way.
            Evan has verified this for beta minus.

        Args:
            beta (str): The beta decay type. Options are ['+', '-', 'c'].
            output (list of str, str): Input for hfbthoParser.

        Returns:
            float: Maximum QRPA energy for beta decay.
        """

        if beta not in [u'+',u'-',u'c']:
            raise ValueError(u"Valid args for kwarg beta are '+', '-', or 'c'.")

        if output is None: output = self.output

        soln = hfbthoParser(output)

        # Get eqrpa_max from the fermi energies
        #   dmass_nh = 0.78227     # neutron hydrogen mass difference
        #   MEC2     = 0.510998928 # electron mass
        #   Be       = 0.0?        # Binding energy of captured electron
        hfb_lambda = soln.getLambda()
        if beta == u'-':
            eqrpa_max = hfb_lambda[0] - hfb_lambda[1] + DMASS_NH
        elif beta  == u'+':
            eqrpa_max = hfb_lambda[1] - hfb_lambda[0] - DMASS_NH - 2*MEC2
        elif beta  == u'c':
            #Zd = self.nucleus[1] - 1 # daughter charge
            #Zk = Zd - 0.35 # Effective charge due to screening from other electrons
            #Bk = 1 - np.sqrt(1 - (alpha*Zk)**2) # Hydrogen-like K-shell BE
            eqrpa_max = hfb_lambda[1] - hfb_lambda[0] - DMASS_NH # - Bk

        return eqrpa_max

    #-----------------------------------------------------------------------
    def getEgs(self, beta, output=None):
        """
        Calculate the transition energy for the transition to daughter ground
        state in the HFB approximation. In this approximation, it is the energy
        of the lowest energy 2 quasiparticle state.

        The energy is:
            * Even-Even: min(EQP_n) + min(EQP_p)
            * Prot-Odd : min(EQP_n) - EQP_blocked_p
            * Neut-Odd : min(EQP_p) - EQP_blocked_n
            * Odd-Odd  : - EQP_blocked_p - EQP_blocked_n
        This approximation has been checked using a global calculation of ~4000
        nuclei, comparing Q = Egs + EQRPAmax to Q = BE_HFB(i) - BE_HFB(f) + dMnH.
        Results cluster right on zero, a typical result may differ by as much as 2 MeV.
        Results with magic numbers (magic +/-1 for odd) agree worse, differing by as
        much as 5-6MeV in some cases.

        Note:
            Recall the QRPA energy contains the HFB Routhian and therfore the fermi
            energy lagrange multiplier. Therfore for charge changing transitions it
            really only reflects the change in pairing energy, not the change in
            particle number of Z and N.

            For an even-even parent the change in pairing energy would be the sum
            of the lowest two quasiparticle energies, as above.

            For odd-A nuclei it would be the minimum 1QP transition energy, which
            would involve removing the odd nucleon (reducing the pairing energy) and
            adding another unpaired nucleon (increasing the pairing energy).

            For odd-odd nuclei we could in principle have the odd states pairing up
            (e.g. a neutron particle fills a proton hole), in which case the pairing
            energy would be reduced by the sum of the two blocked QP energies (I THINK).
            However, in the EFA the occupations forbid this transition, so we don't
            consider it. The next lowest would be of the same form as in the odd-A case.

        Args:
            beta (str): The beta decay type. Options are ['+', '-'].
            output (list of str, str): Input for hfbthoParser.

        Returns:
            float: Ground state transition energy.
        """

        if output is None: output = self.output

        soln = hfbthoParser(output)

        pn_data = soln.getQPData()
        if pn_data is None:
            return np.nan

        # Negate the blocked QP energies - This is a trick to get the HFB
        # (residual int = 0) estimate of Egs in the EFA.
        for pn in list(pn_data.values()):
            if pn[u'blo_qp'] is not None:
                bqp_block = pn[u'blo_qp'][u'block']
                block_ind = pn[u'blo_qp'][u'state']
                bqp_dict  = [d for d in pn[u'qp_dicts'] if d[u'block'] == bqp_block][block_ind]
                for d in pn[u'qp_dicts']:
                    if d == bqp_dict: d[u'eqp'] = -abs(d[u'eqp'])

        # Only consider appropriately occupied/unoccupied qps for beta +/-
        # if pairing has collapsed. Otherwise qp's have partial occupations
        # and transitions are allowed. Note blocked qps have occ=0.5 b/c EFA.
        occ_tol = 5e-9
        if beta == u'-':
            pfac = 1.0; nfac = 0.0
        else:
            pfac = 0.0; nfac = 1.0
        pn_data[u'p'][u'occ_fac'] = pfac
        pn_data[u'n'][u'occ_fac'] = nfac
        for pn in list(pn_data.values()):
            pn[u'qp_dicts'] = [d for d in pn[u'qp_dicts'] if \
                    not abs(pn[u'occ_fac'] - d[u'occ']) < occ_tol ]

        # If there are no QP's, the Egs is undefined... This can happen in
        # finite temperature if the temp is too high.
        if not pn_data[u'p'][u'qp_dicts'] or not pn_data[u'n'][u'qp_dicts']:
            return np.nan

        # Gather 2qp transitions. This is a dictionary b/c
        # historically I looked at the quantum numbers of the final state as well.
        trans = [{u'E2qp':dp[u'eqp'] + dn[u'eqp']}\
                    for dp in pn_data[u'p'][u'qp_dicts']\
                    for dn in pn_data[u'n'][u'qp_dicts']]
                    #if not (dp[u'eqp']<0 and dn[u'eqp']<0)]

        daughter = min(trans, key = lambda x:x[u'E2qp'])

        return daughter[u'E2qp']

    #-----------------------------------------------------------------------
    def getQval(self, beta, **kwargs):
        """
        Calculate the HFB estimate of the beta decay Q-value.

        Value is determined from the approximation in (Engel et al.,
        Phys. Rev. C 60, 014302 (1999)). See getEqrpaMax for details.

        Args:
            beta (str): The beta decay type ['-','+'].

        Returns:
            float: The beta decay Q-value.
        """

        emax = self.getEqrpaMax(beta, **kwargs)
        egs = self.getEgs(beta, **kwargs)
        return emax - egs

    #-----------------------------------------------------------------------
    def getSoln(self, output=None):
        """
        Extract and/or construct key details for the HFB solution and populate
        the soln_dict.

        Note if hfbthoRun.beta is None, values which depend on this will not be
        returned in the soln_dict (Q-value, EQRPA_max, E_gs).
        """

        if output is None: output = self.output

        soln = hfbthoParser(output)
        sd = {}

        # Values parsed directly from output
        energy = soln.getEnergy()
        time   = soln.getTime()
        conv   = soln.getConv()
        defo   = soln.getBeta2()[2]
        sd.update({u'Def':defo,u'Energy':energy,u'Time':time,u'Conv':conv})

        # Values constructed from output (soln.output to avoid reading file again)
        kp = self.getSpinPar(soln.output)
        sd.update({u'KP':kp})

        # Blocking quantum numbers (could be different than user requested)
        blo_out = self.getBlockingResult(soln.output)
        blo_str = self.get_blocking_str(blo_out)
        sd.update({u'Z-Blk': blo_str[1], u'N-Blk': blo_str[0]})

        # Values constructed from output with beta_type dependency
        if self.beta is not None:
            # Q-value related quantities come from 0T !GROUND STATE!. If no such soln can be
            # found, the parser will place nan's, which is my desired result.
            if self.ft_active:
                soln.updateOutput(os.path.join(self.paths.hfb, u'zeroT', self.file_txt))
            emax = self.getEqrpaMax(self.beta, soln.output)
            egs  = self.getEgs(self.beta, soln.output)
            Q = emax - egs # Avoid re-running the parser...
            sd.update({u'HFB_Qval':Q, u'EQRPA_max':emax, u'E_gs':egs})

        # Check this last to yield the correct error message
        soln.findErrors()

        # Store results (should this return instead?)
        return sd, soln.err, soln.err_msg

    #-----------------------------------------------------------------------
    def runExe(self, stdout_bool=False, debug=0):
        """
        Extend the runExe method to handle the hfbtho output.

        Args:
            stdout_bool (bool): Print to stdout in real time.
            debug (int): If !=0 don't actually run the executable.

        Returns:
            str, None : Stderr or None if empty.
        """

        # Returns stdout, stderr.
        stdout, err = fortProcess.runExe(self, stdout_bool=stdout_bool, debug=debug)


        # For Q-val, we need the QP info from thoout.dat not included in stdout
        # So just read from file instead (stdout is thus unused).
        outfile = os.path.join(self.rundir, hfbthoRun.file_txt)
        if not os.path.exists(outfile) and err is None:
            err = u'Missing '+hfbthoRun.file_txt

        # Gather the solution details
        self.updateOutput(outfile, get_soln=True)

        # Write log and collect errors.
        log = hfbLogger(hfbthoRun.file_log)
        log.quickWrite(self.soln_dict, self.rundir)
        if err is None and self.soln_err:
            err = self.soln_err_msg

        # Check errors printed to stdout
        soln = hfbthoParser(stdout)
        soln.findErrors()
        if err is None and soln.err:
            err = soln.err_msg

        # For now, treat blocking cand not found as non-conv
        if soln.err_msg == soln.keys[u'cand_not_found']:
            err = None

        return err
