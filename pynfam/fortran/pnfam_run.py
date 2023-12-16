# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import str
# -------------- Utilities -----------------
from collections import OrderedDict as OD
import numpy as np
import copy
import os
# ----------- Relative Imports -------------
from ..outputs.pnfam_parser import pnfamParser
from .fortran_utils         import fortProcess
from .hfbtho_run            import hfbthoRun
from ..config               import DEFAULTS
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS pnfamRun                                    #
#===============================================================================#
class pnfamRun(fortProcess):
    """
    A pnfamRun is a fortProcess subclass that contains the methods to run the
    serial pnfam fortran program for a single energy point and handle the output.

    All of the higher level beta decay calculations are handled by python. The
    pnfamRun class is effectively one iterative solution of the pnfam equations
    for a single energy point and operator, analogous to one run of HFBTHO.

    Args:
        paths (pynfamPaths): The paths for the pynfam calculation.
        operator (str): The operator name, as in the pnfam namelist.
        kval (int): The operator K-value, as in the pnfam namelist.

    Attributes:
        op (str): The operator name, as in the pnfam namelist.
        k (int): The operator K-value, as in the pnfam namelist.
        soln_dict (dict): Contains key outputs.
        soln_file (file): The pnfam output file name OP.out.
        _soln (pnfamOutParser): Internal initialization for soln attribute.
        _log (pnfamLogParser): Internal initialization for log attribute.
        ctr_params (dict): Contour parameters for the EQRPA point.
    """

    file_exe=u"pnfam_main.x"
    """ file_exe (str): File name of executable.
    """

    def __init__(self, paths, operator, kval):
        # Operator characteristics
        self.op = operator
        self.k = kval

        # Inheret/override attributes for fortran object
        exe = pnfamRun.file_exe
        default_nml = self._getDefaultNml(DEFAULTS[u'fam'], paths)
        fortProcess.__init__(self, paths, exe, self.opname+u'.in', default_nml)

        # Strength result
        self.output = None
        self.soln_err = False
        self.soln_err_msg = u''
        self.soln_dict = {}

        # Contour parameters (Note EQRPA point comes from namelist-see property omega)
        self.ctr_params = None

    # File names based off of operator
    @property
    def file_nml(self):
        """ file_nml (str): File name of namelist.
        """
        return u"{:}.in".format(self.opname)
    @property
    def file_txt(self):
        """ file_txt (str): File name of textfile output.
        """
        return u"{:}.dat".format(self.opname)
    @property
    def file_bin(self):
        """ file_bin (str): File name of binary output.
        """
        return u"{:}.bin".format(self.opname)

    # Different info/strings for the operator
    @property
    def opname(self):
        """ opname (str): Bare_operator + beta_type + kval.
        """
        return u"{:}K{:1d}".format(self.op, self.k)
    @property
    def bareop(self):
        """ bareop (str): Bare operator name without beta type.
        """
        return self.op[:-1]
    @property
    def beta(self):
        """ beta (str): Beta decay type.
        """
        return self.op[-1]
    @property
    def genopname(self):
        """ genopname (str): Bare_operator + _ + kval.
        """
        return u"{:}_K{:1d}".format(self.bareop, self.k)

    # Meta data
    @property
    def omega(self):
        """ omega (complex): The complex energy point.
        """
        return complex(self.nml[u'GENERAL'][u'real_eqrpa'],
                       self.nml[u'GENERAL'][u'imag_eqrpa'])
    @property
    def interaction(self):
        """ interaction (str): The interaction name as in the pnfam namelist.
        """
        return self.nml[u'INTERACTION'][u'interaction_name']



    #-----------------------------------------------------------------------
    def setCtrPoint(self, ctr_pt):
        """
        Translate a complex energy into the appropriate namelist inputs.

        Args:
            ctr_pt (complex): The complex energy point.
        """

        setpoint = {u'GENERAL':
                       {u'real_eqrpa': np.real(ctr_pt),
                        u'imag_eqrpa': np.imag(ctr_pt)}
                    }
        self.adjustNml(setpoint)

    #-----------------------------------------------------------------------
    def runExe(self, stdout_bool=False, debug=0):
        """
        Extend the runExe method to handle the pnfam output.

        Args:
            stdout_bool (bool): Print to stdout in real time.
            debug (int): If !=0 don't actually run the executable.

        Returns:
            str, None : Stderr or None if empty.
        """

        # Returns stdout, stderr.
        stmt = [os.path.join(self.paths.exe, self.exe), self.fname_nml]
        out, err = fortProcess.runExe(self, stdout_bool=stdout_bool,
            statement=stmt, debug=debug)

        outfile = os.path.join(self.rundir, self.file_txt)
        if not os.path.exists(outfile) and err is None:
            err = u"Missing pnFAM output file."

        self.updateOutput(outfile, get_soln=True)

        # Collect errors (output errors are printed over log errors)
        if self.soln_err and err is None:
            err = self.soln_err_msg

        return err

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
    def getSoln(self, output=None):
        """ Parse log file and return soln_dict.
        """

        if output is None: output = self.output

        soln = pnfamParser(output)
        sd = {}

        sd[u'Time'] = soln.getTime()
        sd[u'Conv'] = soln.getConv()
        sd[u'Version']  = soln.getVersion()
        sd[u'Strength'] = soln.getStrength()
        sd[u'Xterms']   = soln.getXterms()

        return sd, soln.err, soln.err_msg

    #-----------------------------------------------------------------------
    def _getDefaultNml(self, DEFAULTS_FAM, paths):
        """
        Combine the default fam settings with the rest of the parameters
        needed to run pnfam but are handled internally.

        Args:
            DEFAULTS_FAM (OrderedDict): The default values from config.py.

        Returns:
            OrderedDict: The full namelist for pnfam.
        """

        # Default dictionary from config.py
        default_nml = copy.deepcopy(DEFAULTS_FAM)

        # Some namelist parameters are not available to the user and are defined
        # interanlly here, rather than in config.py. Some are determined by the
        # pynfam mode (e.g. operators, beta type) while others are set later by
        # other objects, e.g. contour sets the energies.
        fam_internal_nml = OD([
        (u'GENERAL', OD([
            (u'fam_output_filename' , self.opname),
            (u'print_stdout'        , True),
            (u'use_fam_storage'     , None),
            (u'real_eqrpa'          , 0.0),
            (u'imag_eqrpa'          , 0.1)
            ])),
        (u'EXT_FIELD', OD([
            (u'beta_type'           , self.beta),
            (u'operator_name'       , self.bareop),
            (u'operator_k'          , self.k),
            (u'compute_crossterms'  , None)
            ])),
        (u'INTERACTION', OD([
            (u'interaction_name'    , None),
            (u'require_self_consistency', None),
            (u'require_gauge_invariance', None),
            (u'force_j2_terms'      , None),
            (u'vpair_t0'            , None),
            (u'vpair_t1'            , None),
            (u'override_cs0'        , None),
            (u'override_csr'        , None),
            (u'override_cds'        , None),
            (u'override_ct'         , None),
            (u'override_cgs'        , None),
            (u'override_cf'         , None)
            ])),
        (u'SOLVER', OD([
            (u'max_iter'            , None),
            (u'convergence_epsilon' , None),
            (u'broyden_history_size', None),
            (u'energy_shift_prot'   , None),
            (u'energy_shift_neut'   , None),
            (u'quench_residual_int' , None)
            ]))
        ])

        # Override all the None's in the internal namelist with values from config.py
        # (assumes config.py nml has same headers as internal nml)
        for h in fam_internal_nml:
            fam_internal_nml[h].update(default_nml[h])

        return fam_internal_nml
