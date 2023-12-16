# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
# --------------- Utilties -----------------
from shutil import copy2
from shutil import rmtree
from copy import deepcopy, copy
import os
import tarfile
import numpy as np
# ----------- Relative Imports -------------
from .utilities.mpi_utils import *
from .strength.fam_strength import famStrength
from .strength.phase_space import phaseSpace
from .strength.contour import famContour
from .fortran.pnfam_run import pnfamRun
from .utilities.contour_utils import multiple_psi_or_ctr, get_contour_set, contourSetElem
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-08-28'

#-------------------------------------------------------------------------------
def initialize_fam_2bc(mgr, fam_ops, setts, rerun, stdout, nshells):

    nml_params = setts[u'fam']
    tbc_mode = nml_params.get(u'two_body_current_mode', 0)
    if not tbc_mode:
        return
    else:
        tbc_mode = int(str(tbc_mode)[1]) # 2nd digit = 1 needs .tbc
        if tbc_mode != 1:
            return

    err = []
    for d in fam_ops:
        fam_rep = pnfamRun(mgr.paths, d[u'op'], d[u'k'])
        fam_rep.label = mgr.paths.hfb
        if nml_params is not None:
            fam_rep.setNmlParam(nml_params)
        fam_rep.setCtrPoint(1+1j) # Dummy energy
        tbc_input = os.path.join(mgr.paths.hfb, u'{:}.tbc'.format(fam_rep.opname))

        if not os.path.exists(tbc_input):
            if rerun != 2:
                if nshells <= 6: break
                err = [u"Required .tbc file not found in hfb_soln."]
                break
            else:
                err.append(fam_rep.runExe(stdout))

    #**** Check for Errors ****#
    err = [e for e in err if e is not None]
    err_run = bool(err)
    if err_run:
        msg = ["Error encountered initializing 2BC."]+err
        pynfam_warn(msg, mgr.paths.calclabel, err_run)
    #**************************#

    return err_run

#-------------------------------------------------------------------------------
def initialize_fam_contour(mgr, setts, ctr_type, beta_type, hfb_gs):

    # Instantiate and apply overrides
    ctr_main = famContour(ctr_type, setts[u'ctr'])

    # Change the INTERVAL according to hfb_gs (default pynfam behavior)
    shift = 0.0
    for suffix in [u'_prot', u'_neut']:
        try:
            shift += setts[u'fam'][u'energy_shift'+suffix]
        except KeyError:
            pass
    ctr_main.setHfbInterval(hfb_gs, beta_type, shift=shift)

    # Override the INTERVAL if supplied, otherwise this is redundant/has no effect
    ctr_main.updateSettings(setts[u'ctr'])

    ctre = contourSetElem(ctr_main, None, u"")
    contours = [ctre]
    if ctre.ctr_valid:
        # Try since we verify user inputs for psi are valid here
        try:
            contours = multiple_psi_or_ctr(mgr, setts, hfb_gs, ctr_main)
        except RuntimeError as msg:
            pynfam_warn(msg, mgr.paths.calclabel)
            return None

    if all(c.ctr_id is None for c in contours):
        #**** Check for Errors ****#
        msg = [u"Invalid EQRPA interval. Skipping FAM calculation.",
               u"* If using HFB auto-interval, this means Q<=0. To calculate",
               u"  strength anyway, manually specify a valid interval."]
        pynfam_warn(msg, mgr.paths.calclabel)
        return None # FAM errors do not break the dripline
        #**************************#

    return contours

#-------------------------------------------------------------------------------
def initialize_fam_calc(mgr, setts, fam_ops, contours, hfb_gs):
    """
    Populate lists of finished and unfinished fam tasks, and a dictionary
    of famStrength objects for every operator we are computing.
    - Multiple PSI or Contours:
      all_strs_mc has the same form as contours, but famContour objects are replaced
      with dictionaries of famStrength objects on that contour.
    """

    all_strs_set, all_fam_unf, all_fam_fin = [], [], []
    fam_unf, fam_fin = [], []

    contour_set = get_contour_set(contours)
    for ctre in contour_set:
        ctr_main = ctre.ctr

        # Make any subdirectories for multiple (unique) ctrs
        subdir = u""
        if len(contour_set) != 1:
            subdir = ctre.subdir
            ctre.mkdirs(mgr, fam=True)

        all_strengths = {}
        for d in fam_ops:
            strength = famStrength(d[u'op'], d[u'k'], ctr_main, hfb_gs)

            # Make the top level op dir here to avoid possible race condition on calls to runExe?
            mgr.paths.mkdirs(os.path.join(mgr.paths.fam_m, subdir, strength.opname))

            fam_unf, fam_fin = mgr.getFamState(strength, setts[u'fam'], subdir=subdir)

            all_strengths[d[u'name']] = strength
            all_fam_fin += fam_fin
            all_fam_unf += fam_unf

        ctre.all_strengths = all_strengths

    return contour_set, (all_fam_unf, all_fam_fin)

#-------------------------------------------------------------------------------
def run_fam_calc(comm, mgr, stdout, all_fam_unf, all_fam_fin):

    finished_pts, err_run, err_msg = runtasks_master(all_fam_unf, comm, stdout)

    #**** Check for Errors ****#
    msg = [pnfamRun.file_exe+u" encountered an error.",err_msg]
    if pynfam_warn(msg, mgr.paths.calclabel, err_run):
        return None # FAM errors do not break the dripline
    #**************************#

    all_fam_fin += finished_pts
    # Copy a log with EDF info to main dir
    if all_fam_fin:
        o = all_fam_fin[0]
        copy2(os.path.join(o.rundir,o.file_txt),
              os.path.join(mgr.paths.fam_m,u"fam.log"))

    return all_fam_fin

#-------------------------------------------------------------------------------
def finalize_fam_calc(mgr, contour_set, all_fam_fin):

    for ctre in contour_set:
        ctr_main = ctre.ctr
        all_strengths = ctre.all_strengths
        subdir = u""
        if len(contour_set) != 1:
            subdir = ctre.subdir

        for name, strength in list(all_strengths.items()):
            if strength.str_df is None:
                fam_by_op = [o for o in all_fam_fin if o.opname==name and subdir in o.label]

                # We may have lost the order of the points on gather, so sort by label
                fam_by_op.sort(key=lambda x: x.label)

                # Extract strength dataframes from lists of fam objects and write output
                dest = os.path.join(mgr.paths.fam, subdir)
                strength.concatFamData(fam_by_op)
                strength.writeStrengthOut(dest=dest)
                strength.writeCtrBinary(dest=dest)

        # Only tar if we have untarred data, as indicated by all_fam_fin not empty
        # (tar in parallel here, since there could be many directories)
        if all_fam_fin:
            tar_itmp = 0 # newrank
            tar_tasks= [[o for o in all_fam_fin if o.opname==name and subdir in o.label] \
                for name in all_strengths]
            while True:
                if tar_itmp > len(tar_tasks)-1: break
                tar_index = tar_itmp
                tar_itmp += 1 # newcomm_size
                tar_fam   = tar_tasks[tar_index]
                tar_fam_solns(mgr.paths, tar_fam)

# ------------------------------------------------------------------------------
def tar_fam_solns(all_paths, obj_list):
    """
    Bundle individual fam points in a tarfile for portability.

    Args:
        all_paths (pynfamPaths): Paths of the pynfam calculation.
        obj_list (list of pnfamRun): Fortran program run directories to be tarred.
    """
    if not obj_list: return
    top_dir1 = os.path.dirname(obj_list[0].rundir)
    top_dir2 = os.path.dirname(top_dir1)
    tar_name = top_dir1+u'.tar'
    with tarfile.open(tar_name,u'w') as tar:
        for obj in obj_list:
            tar.add(obj.rundir, arcname=os.path.relpath(obj.rundir, top_dir2))
            rmtree(obj.rundir)
    rmtree(top_dir1)
