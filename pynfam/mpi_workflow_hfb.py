# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
import sys
# ----------- Relative Imports -------------
from .utilities.mpi_utils import *
from .utilities.hfb_utils import *
from .fortran.hfbtho_run import hfbthoRun
from .outputs.ls_logger import hfbLogger
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-08-28'


#-------------------------------------------------------------------------------
def initialize_hfb_calc(mgr, beta_type, gs_def_scan, setts, dripline, index2):
    """
    Initialize the main hfbthoRun object template with desired settings, evaluate the
    state of existing outputs, and return finished/unfinished tasks (on all processes).

    Args:
        mgr (pynfamManager): For interfacing with existing outputs and paths.
        beta_type (str): The type of beta decay.
        gs_def_scan (tuple): The deformation parameters from the input file.
        setts (dict): The dictionary of settings for the pynfam calc.
        dripline (int): The dripline mode, as in the input file.
        index2 (int): The dripline loop index.

    Returns:
        hfbthoRun: The main hfbthoRun object template.
        hfbthoRun, None: The ground state, if available from existing ouptuts.
        tuple of hfbthoRun: (Unfinished tasks, finished tasks) for even nuclei.
        tuple of hfbthoRun: (Unfinished tasks, finished tasks) for even nuclei.
    """

    hfb_main = hfbthoRun(mgr.paths, beta_type=beta_type)
    hfb_main.setNmlParam(setts[u'hfb'])
    hfb_main.increment_nucleus(dripline, index2, direction=u'N')

    hfb_gs, even_unf_fin, odd_unf_fin = None, (), ()
    hfb_gs, even_unf_fin, odd_unf_fin= mgr.getHfbState(hfb_main, gs_def_scan, beta_type)

    return hfb_main, hfb_gs, even_unf_fin, odd_unf_fin

#-------------------------------------------------------------------------------
def run_hfb_even_calc(Comm, mgr, ignore_nc, stdout, even_unf, even_fin):

    finished_nucs, err_run, err_msg = runtasks_master(even_unf, Comm, stdout)

    #**** Check for Errors ****#
    msg = [hfbthoRun.file_exe+u" encountered an error running initial solutions.",err_msg]
    if pynfam_warn(msg, mgr.paths.calclabel, err_run):
        return None # break dripline loop, end calc
    #**************************#

    # Retry non-conv with a basis deformation determined by non-conv soln
    rerun_tasks, finished_nucs = hfbEvenRerunList(finished_nucs)
    finished_reruns, err_rerun, err_msg = runtasks_master(rerun_tasks, Comm, stdout)
    finished_nucs += finished_reruns

    #**** Check for Errors ****#
    msg = [hfbthoRun.file_exe+u" encountered an error re-running initial solutions.",err_msg]
    if pynfam_warn(msg, mgr.paths.calclabel, err_rerun):
        return None # break dripline loop, end calc
    #**************************#

    even_fin += finished_nucs
    err_conv = False
    err_conv = hfbConvCheck(even_fin, ignore_nc)[1]
    # Write logfile.dat for all hfb runs before quiting if err
    if err_conv:
        hfb_df = mgr.gatherHfbLogs()
        log = hfbLogger()
        log.quickWrite(hfb_df, dest=mgr.paths.hfb_m)

    #**** Check for Errors ****#
    msg = u"An HFB even solution did not converge."
    if pynfam_warn(msg, mgr.paths.calclabel, err_conv):
        return None # break dripline loop, end calc
    #**************************#

    return even_fin

#-------------------------------------------------------------------------------
def run_hfb_odd_calc(comm, mgr, ignore_nc, stdout, gs_def_scan, hfb_main, even_fin, odd_unf, odd_fin):

    # No secondary hfb calcs
    if not (hfb_main.blocking[0][0] or hfb_main.blocking[1][0]) and not hfb_main.ft_active:
        return odd_fin

    # (odd_fin, odd_unf)=([], []) unless even_fin (existing data) is complete and
    # even_unf=[], in which case one or both odd lists will be filled with something.
    # Starting from scratch, we need to populate the odd_unf task list, treating
    # non-converged even solutions as per the ignore_nonconv input.
    if not odd_unf and not odd_fin:
        if hfb_main.ft_active:
            odd_unf = hfbFTList(even_fin, hfb_main)
        else:
            even_conv = hfbConvCheck(even_fin, ignore_nc)[0]
            odd_unf = hfbOddList(even_conv, hfb_main, gs_def_scan[0])

    finished_nucs, err_run, err_msg = runtasks_master(odd_unf, comm, stdout)

    #**** Check for Errors ****#
    msg = [hfbthoRun.file_exe+u" encountered an error running odd solutions.",err_msg]
    if pynfam_warn(msg, mgr.paths.calclabel, err_run):
        return None # break dripline loop, end calc
    #**************************#

    odd_fin += finished_nucs
    err_conv = False
    err_conv = hfbConvCheck(odd_fin, ignore_nc)[1]
    # Write logfile.dat for all hfb runs before quiting if err
    if err_conv:
        hfb_df = mgr.gatherHfbLogs()
        log = hfbLogger()
        log.quickWrite(hfb_df, dest=mgr.paths.hfb_m)

    #**** Check for Errors ****#
    msg = u"An HFB odd solution did not converge."
    if pynfam_warn(msg, mgr.paths.calclabel, err_conv):
        return None # break dripline loop, end calc
    #**************************#

    return odd_fin

#-------------------------------------------------------------------------------
def finalize_hfb_gs(mgr, ignore_nc, hfb_main, odd_fin, even_fin):

    hfb_gs = None
    # hfb_gs might still be None here if no viable solutions...
    if (hfb_main.blocking[0][0] or hfb_main.blocking[1][0]) or hfb_main.ft_active:
        hfb_fin = odd_fin
    else:
        hfb_fin = even_fin
    hfb_conv = hfbConvCheck(hfb_fin, ignore_nc)[0]
    hfb_gs = hfbGroundState(hfb_conv)

    gs_found = (hfb_gs is not None)

    # Copy the ground state to the hfb solution dir. Also update the output to
    # the new location since we will tar the originals, causing hfb parser to fail
    if gs_found:
        if hfb_main.ft_active != hfb_gs.ft_active:
            zt_dir = u'zeroT'
            dest = os.path.join(mgr.paths.hfb, zt_dir)
            mgr.paths.mkdirs(dest)
        else:
            dest = mgr.paths.hfb
        mgr.paths.copyAllFiles(hfb_gs.rundir, dest)
        hfb_gs.label = mgr.paths.rp(dest)
        new_outfile = os.path.join(dest, hfbthoRun.file_txt)
        hfb_gs.updateOutput(new_outfile, get_soln=False)

    #**** Check for Errors ****#
    msg = u"No solutions were viable ground states."
    if pynfam_warn(msg, mgr.paths.calclabel, (not gs_found)):
        return None # break dripline loop, end calc
    #**************************#

    return hfb_gs

#-------------------------------------------------------------------------------
def finalize_hfb_calc(mgr, hfb_gs, odd_fin, even_fin):

    # Only do this if we have untarred data, as indicated by odd_fin/even_fin not empty
    if (odd_fin+even_fin):
        # Make a log for all hfb runs
        hfb_df = mgr.gatherHfbLogs()
        log = hfbLogger()
        log.quickWrite(hfb_df, dest=mgr.paths.hfb_m)

        if hfb_gs is not None:
            # Replace time with total time in gs mini log
            hfb_gs.soln_dict[u'Total_Time'] = np.sum(hfb_df[u'Time'].values)
            log = hfbLogger(hfbthoRun.file_log)
            log.quickWrite(hfb_gs.soln_dict, hfb_gs.rundir)
            # Tar the hfb meta files for portability. Skip if we didn't get a viable gs
            tar_hfb_solns(mgr.paths, odd_fin+even_fin)

# ------------------------------------------------------------------------------
def tar_hfb_solns(all_paths, obj_list):
    """
    Bundle unused HFBTHO solutions in a tarfile for portability.

    Args:
        all_paths (pynfamPaths): Paths of the pynfam calculation.
        obj_list (list of hfbthoRun): Fortran program run directories to be tarred.
    """
    if not obj_list: return
    # Create tarfle
    tar_name = os.path.join(all_paths.hfb_m, hfbthoRun.file_tar)
    with tarfile.open(tar_name,u'w') as tar:
        for n in obj_list:
            tar.add(n.rundir, arcname=os.path.relpath(n.rundir, all_paths.hfb_m))

    # Clean directories
    for n in obj_list:
        up1dir = os.path.dirname(n.rundir)
        if up1dir == all_paths.hfb_m:
            # even calcs rundir = hfb_meta/xxxx/
            rmtree(n.rundir)
        else:
            # odd calcs rundir = hfb_meta/xxxx/xxxx so remove up1dir
            try:
                rmtree(up1dir)
            except OSError:
                pass
