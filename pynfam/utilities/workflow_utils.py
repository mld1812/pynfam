# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from future.builtins   import str
from future.builtins   import range
# -------------- Utilities -----------------
from shutil import copy2, rmtree
import os
import inspect
import tarfile
import pkg_resources
import numpy as np
# ----------- Relative Imports -------------
from ..fortran.hfbtho_run   import hfbthoRun
from ..fortran.pnfam_run    import pnfamRun
from ..outputs.pynfam_paths import pynfamPaths
from ..pynfam_manager       import pynfamManager
from ..config               import DEFAULTS
from .mpi_utils             import pynfam_mpi_traits
from .mpi_utils             import pynfam_warn
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

# ------------------------------------------------------------------------------
def pynfam_input_check(pynfam_inputs):
    """
    Check user provided pynfam_inputs and generate error messages. Also formats
    some inputs e.g. ctr_type --> ctr_type.upper() or None.

    Args:
        pynfam_inputs (dict): The variable from the input script.

    Returns:
        list of str: Error messages.
    """

    # Error messages
    err = []

    # Unpack dicts for convenience
    dirs          = pynfam_inputs[u'directories']
    nr_masters    = pynfam_inputs[u'nr_parallel_calcs']
    rerun         = pynfam_inputs[u'rerun_mode']

    dripline      = pynfam_inputs[u'hfb_mode'][u'dripline_mode']
    gs_def_scan   = pynfam_inputs[u'hfb_mode'][u'gs_def_scan']
    ignore_nc     = pynfam_inputs[u'hfb_mode'][u'ignore_nonconv']

    beta_type     = pynfam_inputs[u'fam_mode'][u'beta_type']
    fam_ops_in    = pynfam_inputs[u'fam_mode'][u'fam_ops']
    ctr_type      = pynfam_inputs[u'fam_mode'][u'fam_contour']

    # Format the ctr_type
    if ctr_type is not None:
        ctr_type = ctr_type.upper()
        if ctr_type in [u'NONE',u'']: ctr_type = None
    pynfam_inputs[u'fam_mode'][u'fam_contour'] = ctr_type

    # Check inputs
    if any(not isinstance(d,str) for d in list(dirs.values())):
        err.append(u"Invalid input for directories. Must be text.")
    if nr_masters is not None and not (isinstance(nr_masters, int) and nr_masters>=0):
        err.append(u"Invalid input for nr_parallel_calcs. Must be positive integer or None")
    if rerun not in [0,1,2]:
        err.append(u"Invalid input for rerun_mode. Options are [0,1,2].")
    if abs(dripline) not in [0,1,2]:
        err.append(u"Invalid input for dripline_mode. Options are [0,+/-1,+/-2].")
    if ignore_nc not in [0,1,2,3]:
        err.append(u"Invalid input for ignore_nonconv. Options are [0,1,2,3].")
    if beta_type not in [u'+', u'-', u'c']:
        err.append(u"Invalid input for beta_type. Options are ['+','-','c'].")
    if not (isinstance(gs_def_scan, tuple) and len(gs_def_scan)==2):
        err.append(u"Invalid form for gs_def_scan. Expected tuple of length 2.")
    else:
        if not isinstance(gs_def_scan[0],int) or not isinstance(gs_def_scan[1],tuple):
            err.append(u"Invalid types for gs_def_scan elements. Expected (int, tuple).")
        elif abs(gs_def_scan[0]) not in [0,1,2]:
            err.append(u"Invalid input for gs_def_scan[0]. Options are [0,+/-1,+/-2].")
    if not isinstance(fam_ops_in, (tuple,str)):
        err.append(u"Invalid input for fam_ops. See documentation.")

    # Check contour type
    ctr_options = [k for k in list(DEFAULTS[u'ctr'].keys()) if k != u'INTERVAL']+[None]
    if ctr_type not in ctr_options:
        err.append(u"Invalid input for ctr_type. See documentation for options.")

    # Format def scan (if kickoff != 0, but tuple is empty, turn off scan)
    if not gs_def_scan[1]:
        pynfam_inputs[u'hfb_mode'][u'gs_def_scan'] = (0, ())


    return err

# ------------------------------------------------------------------------------
def pynfam_override_check(pynfam_inputs, override_settings):
    """
    Check user provided override_settings and generate error messages.

    Args:
        pynfam_inputs (dict): The variable from the input script.
        override_settings (dict): The variable from the input script.

    Returns:
        list of str: Error messages.
    """
    # NB!!! - All override settings can be tuples, so check accordingly!

    f_conf = "config.py"

    # Error messages
    err = []

    # Unpack
    ctr_type = pynfam_inputs[u'fam_mode'][u'fam_contour']

    ctr_set = override_settings[u'ctr']
    psi_set = override_settings[u'psi']

    # Contour settings
    if ctr_type is not None:
        # Check for valid keys
        ctr_defaults  = list(DEFAULTS[u'ctr'][ctr_type].keys())
        ctr_defaults += list(DEFAULTS[u'ctr'][u'INTERVAL'].keys())
        if any(k not in ctr_defaults for k in list(ctr_set.keys())):
            err_str = u"Invalid override setting for ctr: {:}. Check {:} for options."
            err.append(err_str.format(ctr_type,f_conf))
        # Check bounds are real
        for k in ["energy_min", "energy_max"]:
            if k in list(ctr_set.keys()):
                if not (np.array(np.isreal(ctr_set[k])).all()):
                    err.append(u"Invalid contour bounds. Must be real.")
        if "half_width" in list(ctr_set.keys()):
            hw = np.array(ctr_set["half_width"])
            if ((np.logical_not(np.isreal(hw)))&(abs(hw.real)>1e-10)).any():
                err.append(u"Invalid input for contour setting 'half_width'. Must be real or imaginary.")

    # Check the psi_approx type is one that's programmed
    psi_options = [k for k in list(DEFAULTS[u'psi'].keys()) if k != u'PSI']
    if u'psi_approx' in list(psi_set.keys()):
        psi_approx = np.array(psi_set[u'psi_approx'], copy=False, ndmin=1)
        if any(k not in psi_options for k in psi_approx):
            err_str = u"Invalid override setting for 'psi_approx'. Check {:} for options."
            err.append(err_str.format(f_conf))
    # Check other override settings are programmed
    psi_defaults = [hh for h in DEFAULTS[u'psi'] for hh in DEFAULTS[u'psi'][h]]
    if any(k not in psi_defaults for k in list(psi_set.keys())):
        err_str = u"Invalid override setting provided for psi. Check {:} for options."
        err.append(err_str.format(f_conf))

    # HFBTHO and PNFAM namelist inputs
    for t in [u'hfb',u'fam']:
        t_defaults = [k2 for k1 in list(DEFAULTS[t].keys()) for k2 in list(DEFAULTS[t][k1].keys())]
        if any(k not in t_defaults for k in list(override_settings[t].keys())):
            err_str = u"Invalid override setting for {:}. Check {:} for options."
            err.append(err_str.format(t,f_conf))

    # FT EC mode
    if pynfam_inputs[u'fam_mode'][u'beta_type'] == u'c' and ctr_type is not None:
        try:
            if override_settings[u'hfb'][u'set_temperature']:
                psi_approx = psi_set.get(u'psi_approx', DEFAULTS[u'psi']['PSI']['psi_approx'])
                if psi_approx != u"POLYFIT":
                    err.append(u"FT-EC calculations must use psi_approx='POLYFIT'")
        except KeyError:
            pass

    return err

# ------------------------------------------------------------------------------
def pynfam_init(pynfam_inputs, override_settings):
    """
    Check user provided inputs are compatible with pynfam functionality and
    generate error messages. Also generate values describing the specified
    run mode, including number of calculations, whether FAM is active, and
    whether there are existing outputs to be used.

    Args:
        pynfam_inputs (dict): The variable from the input script.
        override_settings (dict): The variable from the input script.

    Returns:
        list of str: Error messages.
        list of str: Warning messages.
        int: Number of pynfam calculations.
        bool: Perform the FAM portion of the calculation.
        list of str: Labels of existing outputs.
    """

    # Error messages
    err, warn = [], []

    # Unpack
    dirs          = pynfam_inputs[u'directories']
    rerun         = pynfam_inputs[u'rerun_mode']
    dripline      = pynfam_inputs[u'hfb_mode'][u'dripline_mode']

    beta_type     = pynfam_inputs[u'fam_mode'][u'beta_type']
    fam_ops_in    = pynfam_inputs[u'fam_mode'][u'fam_ops']
    ctr_type      = pynfam_inputs[u'fam_mode'][u'fam_contour']

    hfb_set = override_settings[u'hfb']
    fam_set = override_settings[u'fam']
    ctr_set = override_settings[u'ctr']
    psi_set = override_settings[u'psi']

    # Executables
    paths = pynfamPaths(u'', dirs[u'outputs'], dirs[u'exes'], dirs[u'scratch'])
    fam_exe = os.path.join(paths.exe, pnfamRun.file_exe)
    hfb_exe = os.path.join(paths.exe, hfbthoRun.file_exe)
    if rerun == 0:
        if not os.path.exists(fam_exe) or not os.path.exists(hfb_exe):
            err.append(u"A fortran executable was not found. Requested rerun mode requires hfbtho and pnfam.")
    else:
        if not os.path.exists(fam_exe):
            err.append(u"A fortran executable was not found. Requested rerun mode requires pnfam only.")

    # FAM active?
    do_fam = True
    if (ctr_type is None) or (not get_fam_ops(fam_ops_in, beta_type)):
        do_fam = False
        warn.append(u"No FAM contour or operators requested.\n"\
                    u"    Calculation will stop after hfb.")

    # Number of calculations
    nr_calcs, tup = None, True
    for sub in [hfb_set, fam_set, ctr_set, psi_set]:
        for h in sub:
            if isinstance(sub[h], tuple):
                if len(sub[h]) == 0:
                    err.append(u"Input tuple with length zero detected")
                if nr_calcs is None:
                    nr_calcs = len(sub[h])
                else:
                    if len(sub[h]) != nr_calcs:
                        err.append(u"Input tuples must have same length")
    if nr_calcs is None:
        nr_calcs, tup = 1, False

    # Existing outputs: if so and rerun=1, nr_calcs determined by how many
    exst_data = []
    if rerun!=0:
        if os.path.exists(paths.out):
            # Get just the labels (not full path)
            exst_data = paths.getLabels(prefix=False)

        if exst_data:
            if len(exst_data) != nr_calcs and tup:
                err.append(u"Number of existing directories does not match "+\
                           u"length of tuple inputs.\n")
            elif not tup:
                nr_calcs = len(exst_data)

            # If existing dirs, make sure they all have hfb ground states
            hfb_solns = [os.path.join(paths.out, p, paths.subdir_hfb) for p in exst_data]
            missing = [p for p in hfb_solns if not pynfamManager.checkHfbthoFiles(p)]
            if missing:
                err.append(u"All existing data must include all ground state solution files for rerun_mode!=0. Missing:")
                err.append(u"\n".join(missing))

            # Dripline mode doesn't make sense in this case
            if dripline:
                err.append(u"Cannot rerun with rerun_mode!=0 and dripline active.")
        elif rerun == 2:
            err.append(u"Cannot rerun with rerun_mode==2 from scratch.")

    return err, warn, nr_calcs, do_fam, exst_data

# ------------------------------------------------------------------------------
def pynfam_mpi_init(pynfam_inputs, nr_calcs, comm, check):

    err, warn = [], []
    rank, comm_size = pynfam_mpi_traits(comm)

    input_name = u"nr_parallel_calcs"
    nr_masters = pynfam_inputs[input_name]
    rerun      = pynfam_inputs[u'rerun_mode']

    # Set nr_masters here if nr_parallel_calcs is set to use defaults or inconsistent.
    warn_str = None
    if nr_masters:
        if nr_masters > nr_calcs:
            warn_str = u"Supplied {:}={:}, but nr_calcs={:}. Using {:}={:}."
            warn_str = warn_str.format(input_name, nr_masters, nr_calcs, input_name, nr_calcs)
            nr_masters = nr_calcs
        elif comm_size == 1 and nr_masters > 1:
            warn_str = u"Supplied {:}={:} for serial calculation. Using {:}=1."
            warn_str = warn_str.format(input_name, nr_masters, input_name)
            nr_masters = 1
    else:
        # False type (False, None, 0) invokes default value without warning message
        nr_masters = nr_calcs
    if warn_str: warn.append(warn_str)


    # Split into master/worker groups
    if comm_size > 1 and rerun == 2:
        stdout = False
        newcomm = comm
        group = 0

    elif comm_size > 1:
        stdout = False

        if rank in range(1,nr_masters+1):
            group = 0 # Masters, ranks 1 to nr_masters
        else:
            group = 1 # Workers, with world rank 0 as leader

        newcomm = comm.Split(group,1)

        # Check we have enough resources
        nr_workers = comm_size - nr_masters - 1
        if nr_workers <= nr_masters and not check:
            e1=u"Too few mpi process reserved. Must have at least 1 worker per parallel calc + 1 lead worker"
            e2=u"Requested: {:}={:}, comm_size={:} --> available_workers={:}".format(input_name, nr_masters, comm_size, nr_workers)
            err.append(e1)
            err.append(e2)
        # Warn if nr_workers ~ nr_masters
        if nr_workers <= nr_masters*1.5:
            warn_str = u"Number of worker processes ({:}) is close to number of"+\
                       u" master processes ({:}). Consider requesting more resources."
            warn.append(warn_str.format(nr_workers, nr_masters))
    else:
        stdout = True
        group = 0
        newcomm = 1

    if check:
        warn.append(u"nr_calcs = {:}".format(nr_calcs))

    return err, warn, newcomm, group, stdout

# ------------------------------------------------------------------------------
def pynfam_init_dirs(pynfam_inputs):
    """
    Initialize the output directory tree and copy meta data. This must only
    be called by a single master rank (COMM_WORLD.Get_rank()=0).

    Args:
        pynfam_inputs (dict): The variable from the input script.
    """

    # Make highest level output dir and output meta
    dirs = pynfam_inputs[u'directories']
    init_paths = pynfamPaths(u'', dirs[u'outputs'], dirs[u'exes'], dirs[u'scratch'])
    init_paths.mkdirs()

    # Copy input script
    frame  = inspect.stack()[-1] # highest level calling script
    module = inspect.getmodule(frame[0])
    fname  = module.__file__
    copy2(fname, init_paths.meta)

    # Copy default settings
    rsrc_pkg  = u'pynfam'
    rsrc_path = u'config.py'
    conf_file = pkg_resources.resource_stream(rsrc_pkg, rsrc_path)
    copy2(conf_file.name, init_paths.meta)

    # Copy the (empty) slurm file so the output is associated with a job number
    copy_slurm(init_paths)

# ------------------------------------------------------------------------------
def pynfam_finalize(pynfam_inputs):
    """
    Copy meta data and construct master logfile once calculation is finished.

    Args:
        pynfam_inputs (dict): The variable from the input script.
    """
    dirs = pynfam_inputs[u'directories']
    fin_paths = pynfamPaths(u'', dirs[u'outputs'], dirs[u'exes'], dirs[u'scratch'])
    mgr = pynfamManager(fin_paths)

    # Copy populated slurm file
    copy_slurm(fin_paths)

    # Write master log to meta
    mgr.gatherMasterLog()
    mgr.gatherMasterLog(zeroed=True)

# ------------------------------------------------------------------------------
def copy_slurm(paths):
    """
    Copy slurm output to meta data directory.

    Args:
        paths (pynfamPaths): The calculation paths.
    """

    slurm_id = os.getenv(u'SLURM_JOB_ID')
    if slurm_id is not None:
        try:
            slurm_file = u'slurm_'+slurm_id+u'.out'
            copy2(os.path.join(paths.cwd,slurm_file),
                    os.path.join(paths.meta, slurm_file))
        except IOError:
            msg=u"Found $(SLURM_JOB_ID)= "+slurm_id+u", but corresponding\n"+\
                u"output file missing. Skipping copy to meta."
            pynfam_warn(msg)

# ------------------------------------------------------------------------------
def pynfam_settings(override_settings, index):
    """
    Extract settings for the "index"th pynfam run based on tuple overrides.

     Args:
        override_settings (dict): The variable from the input script.
        index (int): The index for the specific pynfam calculation.
    """
    ovr = {}
    for s in override_settings:
        sett = override_settings[s]
        ovr[s] = {}
        for h in sett:
            if isinstance(sett[h], tuple):
                ovr[s][h] = sett[h][index]
            else:
                ovr[s][h] = sett[h]
    return ovr

# ------------------------------------------------------------------------------
def get_fam_ops(key, beta_type):
    r"""
    Translate input file FAM operator keywords into list of operators for
    pnfam fortran program.

    Note:
        The following operators are implemented in pnfam, and are separated by
        how they change the :math:`K=\{-J,\dots,J\}` quantum number.

        * 'F'   : (J=0) Fermi operator (identity)
        * 'GT'  : (J=1) Gamow-Teller operator (pauli matrix sigma)
        * 'R'   : (J=1) position vector (r)
        * 'P'   : (J=1) momentum over i (p/i)
        * 'RS0' : (J=0) r and sigma coupled to J=0
        * 'RS1' : (J=1) r and sigma coupled to J=1
        * 'RS2' : (J=2) r and sigma coupled to J=2
        * 'PS0' : (J=0) p/i dotted with sigma (note: this is not equal to
           p/i and sigma coupled to J=0)

        This package wraps sets of operators for convenient beta decay calculations.
        Available calculations include total, total allowed, total first-forbidden, or
        transition by :math:`J^{\pi}`. Calculation names, number of strength functions
        computed (n), and operators considered, are listed below. (note only
        :math:`K \ge 0` operators are computed):

        * 'All'       : n=14, all operators
        * 'Allowed'   : n= 3, [F,GT]
        * 'Forbidden' : n=11, [R,P,RS0,RS1,RS2,PS0]
        * '0+'        : n= 1, [F]
        * '1+'        : n= 2, [GT]
        * '0-'        : n= 2, [RS0,PS0]
        * '1-'        : n= 6, [R,P,RS1]
        * '2-'        : n= 3, [RS2]
        * (Operator,K): n= 1, manually specify 1 strength function with a tuple of form
          (str, int) where str is one of the 8 operators listed above and int is the K value.

    Args:
        key (str): FAM operators to calculate, as in the input file.
        beta_type (str): Type of beta decay.

    Returns:
        list of dict: Active FAM operators.
    """

    # Implemented Operators - {OP : (J, pi)}
    fam_ops = {u'F'    : (0,+1),
               u'GT'   : (1,+1),
               u'P'    : (1,-1),
               u'PS0'  : (0,-1),
               u'R'    : (1,-1),
               u'RS0'  : (0,-1),
               u'RS1'  : (1,-1),
               u'RS2'  : (2,-1)
               }
    # Wrapper keys - {key : (J's, pi's)}
    keys = {u'ALL'      : ([0,1,2], [+1,-1]),
            u'ALLOWED'  : ([0,1,2], [+1]),
            u'FORBIDDEN': ([0,1,2], [-1]),
            u'0+'       : ([0]    , [+1]),
            u'GAMOWTELLER':([1]   , [+1]),
            u'1+'       : ([1]    , [+1]),
            u'0-'       : ([0]    , [-1]),
            u'1-'       : ([1]    , [-1]),
            u'2-'       : ([2]    , [-1])
            }

    if beta_type == u'c':
        beta_type = u'+'

    if isinstance(key, tuple):
        if not key:
            return []
        op = key[0].upper()
        kval = key[1]

        try:
            j = fam_ops[op][0]
        except KeyError:
            raise ValueError(u"Requested FAM operator not implemeneted. Options: '{:}'.".format(\
                    u"', '".join(list(fam_ops.keys()))))
        if kval not in list(range(-j,j+1)):
            raise ValueError(u"Invalid operator/K-value assignment. See documentation.")

        active_ops = [{u'name' : u'{:}{:}K{:1d}'.format(op,beta_type,kval),
                       u'genop': u'{:}_K{:1d}'.format(op,kval),
                       u'k'    : kval,
                       u'op'   : u'{:}{:}'.format(op,beta_type)}]
    else:
        if key is not None:
            key = key.upper()
        if key in [None, u'NONE', u'']:
            return []

        try:
            active_ops = [{u'name' : u'{:}{:}K{:1d}'.format(op,beta_type,kval),
                           u'genop': u'{:}_K{:1d}'.format(op,kval),
                           u'k'    : kval,
                           u'op'   : u'{:}{:}'.format(op,beta_type)} \
                for op, (j,p) in list(fam_ops.items()) for kval in range(0,j+1)\
                    if j in keys[key][0] and p in keys[key][1]]
        except KeyError:
            raise ValueError(u"Invalid FAM operator key. Options: '{:}', or (Operator, K).".format(\
                    u"', '".join(list(keys.keys()))))

    return active_ops
