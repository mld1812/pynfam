# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import range
from builtins   import str
# ---------------- Utilities ---------------
from shutil import copy2
import numpy as np
import os
import sys
from copy import deepcopy
# ----------- Relative Imports -------------
from .pynfam_manager        import pynfamManager
from .outputs.pynfam_paths  import pynfamPaths
from .outputs.ls_logger     import hfbLogger
from .outputs.ls_logger     import betaLogger
from .fortran.hfbtho_run    import hfbthoRun
from .fortran.pnfam_run     import pnfamRun
from .strength.contour      import famContour
from .strength.fam_strength import famStrength
from .strength.phase_space  import phaseSpace
from .strength.shape_factor import shapeFactor
from .utilities.mpi_utils import do_mpi
from .utilities.hfb_utils import convString
from .utilities import mpi_utils as mu
from .utilities import workflow_utils as wu
from .utilities.contour_utils import get_contour_split_set, contourSetElem
from .mpi_workflow_hfb  import *
from .mpi_workflow_fam  import *
from .mpi_workflow_beta import *
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'


#===============================================================================#
#                       Primary Pynfam Functions                                #
#===============================================================================#
def pynfam_mpi_calc(pynfam_inputs, override_settings, check=False):
    """
    A top level function to handle most of the MPI of a large scale pynfam calculation.

    Args:
        pynfam_inputs (dict)
        override_settings (dict)
        check (bool)
    """

    # Setup MPI variables
    if do_mpi:
        comm = mu.MPI.COMM_WORLD
    else:
        comm = 0
    rank, comm_size = mu.pynfam_mpi_traits(comm)

    # User input values valid? (All processes must run this b/c we format some inputs)
    inp_err_i = wu.pynfam_input_check(pynfam_inputs)
    inp_err_o = wu.pynfam_override_check(pynfam_inputs, override_settings)
    inp_err = inp_err_i + inp_err_o
    if rank == 0:
        if inp_err: mu.pynfam_abort(comm, inp_err)

    # Inputs compatible with requested mode? If so retrieve some run parameters.
    inp_err, inp_wrn, nr_calcs, do_fam, exst_data = wu.pynfam_init(pynfam_inputs, override_settings)
    if rank == 0:
        if inp_err: mu.pynfam_abort(comm, inp_err)

    # Setup MPI related parameters, split comm into master/worker.
    inp_err, mpi_wrn, newcomm, group, stdout = wu.pynfam_mpi_init(pynfam_inputs, nr_calcs, comm, check)
    if rank == 0:
        if inp_err and not check: mu.pynfam_abort(comm, inp_err)
        if mpi_wrn: mu.pynfam_warn(mpi_wrn)

    # STOP here if we just wanted to check all the inputs but not run.
    if rank==0 and check:
        if inp_wrn: mu.pynfam_warn(inp_wrn)
        check_msg = u"Input check complete. No errors were detected."
        mu.pynfam_abort(comm, check_msg)

    # Setup PyNFAM directory tree
    if rank==0:
        wu.pynfam_init_dirs(pynfam_inputs)

    # Assign calcs to masters round-robin
    if group == 0:
        kwargs = {u'pynfam_inputs': pynfam_inputs,
                  u'override_settings': override_settings,
                  u'comm': comm,
                  u'do_fam': do_fam,
                  u'stdout': stdout,
                  u'exst_data': exst_data
                  }
        rank, nr_masters = mu.pynfam_mpi_traits(newcomm)
        icalc = rank
        while True:
            if icalc > nr_calcs - 1: break
            index = icalc
            icalc += nr_masters
            kwargs.update({u'index': index})
            pynfam_dripline_calc(**kwargs)

        if do_mpi: mu.runtasks_killsignal(comm, newcomm) # Barrier
    # All else are workers. They don't know anything, just run tasks
    else:
        mu.runtasks_worker(comm, newcomm, stdout)

    if do_mpi: comm.Barrier()
    if rank==0: wu.pynfam_finalize(pynfam_inputs)


# ------------------------------------------------------------------------------
def pynfam_dripline_calc(pynfam_inputs, **kwargs):
    r"""
    Wrapper for the main pynfam_calc function which loops until the dripline
    has been reached, as indicated by the one or two nucleon separation energy
    becoming negative.

    .. math::
        Sn_x &= B(N-x,Z) - B(N,Z)\\
        Sp_x &= B(N,Z-x) - B(N,Z)

    Args:
        pynfam_inputs (dict): pynfam_inputs as in the main input file.
        **kwargs
    """

    # Initialializations
    i_drip = 0
    tEnergy_prev = None
    dripline = pynfam_inputs[u'hfb_mode'][u'dripline_mode']

    # Dripline Loop
    while True:
        kwargs.update({
            u'pynfam_inputs':pynfam_inputs,
            u'index2': i_drip
            })
        hfb_gs = pynfam_calc(**kwargs)

        if not dripline or hfb_gs is None:
            break
        if tEnergy_prev is not None:
            dE = tEnergy_prev - hfb_gs.soln_dict[u'Energy']
            if dE < 0.0:
                break
        tEnergy_prev = hfb_gs.soln_dict[u'Energy']
        i_drip += 1

# ------------------------------------------------------------------------------
def pynfam_calc(index, index2, pynfam_inputs, override_settings, comm, do_fam, stdout, exst_data):
    """
    The main pynfam calculation.

    The calculation is carried out in 3 stages:
        1) Find the HFB ground state (deformed, even or odd)
        2) Compute the FAM strength along a contour (open or closed contour)
        3) Integrate phase space weighted strength for rates

    Args:
        pynfam_inputs (dict)
        override_settings (dict)
        comm (mpi4py.MPI.intracomm, int)
        do_fam (bool)
        stdout (bool)
        exst_data (list)
    """
    #---------------------------------------------------------------------------#
    #                        Initialize PyNFAM Calc                             #
    #---------------------------------------------------------------------------#
    dirs        = pynfam_inputs[u'directories']
    rerun       = pynfam_inputs[u'rerun_mode']

    dripline    = pynfam_inputs[u'hfb_mode']['dripline_mode']
    gs_def_scan = pynfam_inputs[u'hfb_mode']['gs_def_scan']
    if isinstance(gs_def_scan, list):
        gs_def_scan = gs_def_scan[index]
    ignore_nc   = pynfam_inputs[u'hfb_mode']['ignore_nonconv']

    beta_type   = pynfam_inputs[u'fam_mode']['beta_type']
    if isinstance(beta_type, list):
        beta_type = beta_type[index]
    fam_ops_in  = pynfam_inputs[u'fam_mode']['fam_ops']
    if isinstance(fam_ops_in, list):
        fam_ops_in = fam_ops_in[index]
    ctr_type    = pynfam_inputs[u'fam_mode']['fam_contour']

    fam_ops = wu.get_fam_ops(fam_ops_in, beta_type)
    setts = wu.pynfam_settings(override_settings, index)

    # Unique run label
    if rerun!=0 and exst_data:
        calclabel = exst_data[index]
    else:
        calclabel = str(index).zfill(6)
    if dripline:
        calclabel_base = calclabel
        calclabel += u'_'+str(index2).zfill(4)

    # Paths (requires a calclabel, otherwise not fully defined) and manager
    all_paths = pynfamPaths(calclabel, dirs[u'outputs'], dirs[u'exes'], dirs[u'scratch'])
    mgr = pynfamManager(all_paths)
    mgr.paths.mkdirs()

    #---------------------------------------------------------------------------#
    #                               HFB Calc                                    #
    #---------------------------------------------------------------------------#
    hfb_main, hfb_gs, hfb_evens, hfb_odds = initialize_hfb_calc(mgr,
            beta_type, gs_def_scan, setts, dripline, index2)
    # In case we have hfb_gs already, but need to tar files, set fin lists
    even_fin, odd_fin = hfb_evens[1], hfb_odds[1]

    if hfb_gs is None:
        even_fin = run_hfb_even_calc(comm, mgr, ignore_nc, stdout, *hfb_evens)
        if even_fin is None: return

        # We need the 0T gs for FT calcs logfiles. If no viable gs, end calc
        if hfb_main.ft_active:
            hfb_gs_0T = finalize_hfb_gs(mgr, ignore_nc, hfb_main, odd_fin=even_fin, even_fin=[])
            if hfb_gs_0T is None: return

        odd_fin = run_hfb_odd_calc(comm, mgr, ignore_nc, stdout,
                gs_def_scan, hfb_main, even_fin, *hfb_odds)
        if odd_fin is None: return

        hfb_gs = finalize_hfb_gs(mgr, ignore_nc, hfb_main, odd_fin, even_fin)

    finalize_hfb_calc(mgr, hfb_gs, odd_fin, even_fin)
    if hfb_gs is None: return

    nshells = hfb_gs.nml[u'HFBTHO_GENERAL'][u'number_of_shells']

    err = initialize_fam_2bc(comm, mgr, fam_ops, setts, rerun, stdout, nshells)
    if rerun == 2 or err:
        if not err:
            msg = u"Rerun mode for 2BC calculation requested. Exiting after 2BC calc."
            mu.pynfam_warn(msg, calclabel)
        return # Break dripline, should not be doing hfb calcs for this mode

    if not do_fam:
        msg = u"No FAM contour or operators requested. Exiting after hfb calc."
        mu.pynfam_warn(msg, calclabel)
        return hfb_gs

    #------------------------------------------------------------------------#
    #                               FAM Calc                                 #
    #------------------------------------------------------------------------#
    ctr_main = initialize_fam_contour(mgr, setts, ctr_type, beta_type, hfb_gs)
    if ctr_main is None: return hfb_gs

    contour_set, fam_pts = initialize_fam_calc(mgr, setts, fam_ops, ctr_main, hfb_gs)

    all_fam_fin = run_fam_calc(comm, mgr, stdout, *fam_pts)
    if all_fam_fin is None: return hfb_gs

    finalize_fam_calc(mgr, contour_set, all_fam_fin)

    #------------------------------------------------------------------------#
    #                               Beta Calc                                #
    #------------------------------------------------------------------------#
    check_qval(hfb_gs)

    try:
        beta_logs = []
        all_strs_set = {c.ctr_id: c.all_strengths for c in contour_set}
        ctr_split_set = get_contour_split_set(ctr_main)
        #print(f"Num elements in ctr split set: f{len(ctr_split_set)}") #this should only have 1 element when no psi overrides - confirmed.
        for split_list in ctr_split_set:
            tot_rates_df  = None
            tot_ratesz_df = None
            #print(f"Num elements in split_list: {len(split_list)}") #also only 1 element - confirmed.
            for ctre in split_list: #for one contour, we can ignore the outer two for-loops
                if ctre.ctr_id is None:
                    ctre.mkdirs(mgr, fam=False)
                    continue
                if ctre.is_ft_pole:
                    continue

                ctre.all_strengths = deepcopy(all_strs_set[ctre.ctr_id])
                closed = ctre.ctr.closed

                # Set phase space settings, handling list inputs.
                # (We use the same psi settings for every split contour)
                psi_setts = deepcopy(setts[u'psi'])
                psi_log = {}
                Evloss = None
                # If list, label != "", set the param and add to logfile
                if ctre.label1:
                    for h in setts[u'psi']:
                        if isinstance(setts[u'psi'][h], list):
                            param = setts[u'psi'][h][ctre.index]
                            psi_setts[h] = param
                            psi_log[h] = param
                if u'log(pYe)' in setts[u'psi']: Evloss = hfb_gs


                # Every contour touching zero has an FT_pole subdir (since phase space is unique to contour)
                btot_ft_pole = None
                rates_ft_pole = None
                pole_time = 0.0
                lab3_list = [u'']
                if ctre.touching_ft_pole:
                    lab3_list = [contourSetElem.FT_pole_label] + lab3_list

                for lab3 in lab3_list: #if Finite temp is not active, then we only have lab3 = ''.
                    # Make the beta directories
                    lab = os.path.normpath(os.path.join(ctre.subdir,lab3))
                    ctre.mkdirs(mgr, fam=False, label3=lab3)

                    # Use the correct strength function for pole vs not-pole
                    if lab3:
                        all_strengths_use = all_strs_set[contourSetElem.FT_pole_label]
                    else:
                        all_strengths_use = ctre.all_strengths

                    # Use the not-pole contour interval for phase space approx for pole and not-pole
                    #print(f"Length of strengths: {len(all_strengths_use)}, key: {all_strengths_use.keys()}") #TESTSAVE
                    shapefacs = shapeFactor(list(all_strengths_use.values()), beta_type, ctre.ctr)
                    shapefacs.updateSettings(psi_setts)

                    # Calculate and write raw outputs (direct result from pnfam)
                    shapefacs.calcShapeFactor(hfb_gs, zero_neg=False)
                    btot_df = shapefacs.calcTotalStr(zero_neg=False)
                    rates_df = calc_beta_rates(shapefacs, closed, hfb_gs=Evloss, zero_neg_str=False)

                    # Loop over pole 1st, then subtract it
                    if lab3:
                        btot_ft_pole = btot_df.div(2.0)
                        rates_ft_pole = rates_df.div(2.0)
                        pole_time = shapefacs.total_time/len([c for c in ctr_main if c.touching_ft_pole])
                    elif rates_ft_pole is not None:
                        btot_df = btot_df.subtract(btot_ft_pole)
                        rates_df = rates_df.subtract(rates_ft_pole)
                        rates_df = shapeFactor.recalcHalflives(rates_df)

                    results_beta_rates(mgr, shapefacs, closed, rates_df, zero_neg=False, subdir=lab)
                    results_shape_factor(mgr, shapefacs, closed, write_psi=True, subdir=lab)
                    results_bare_strength(mgr, shapefacs, closed, btot_df=btot_df, subdir=lab)
                    write_beta_log(mgr, shapefacs, rates_df, hfb_gs, psi_log, subdir=lab)

                    # Calculate and write adjusted outputs (written to beta_meta).
                    # 1. Zero negative rate contribs in beta.out
                    # 2. Shapefactor and Btot calculated with zeroed negative strength
                    #    Note: this uses findFirstPeak (which may fail) to set all str
                    #          before 1st peak to zero.
                    if not closed:
                        shapefacs = shapeFactor(list(all_strengths_use.values()), beta_type)
                        shapefacs.updateSettings(psi_setts)
                        shapefacs.calcShapeFactor(hfb_gs, zero_neg=True)

                        rates_df = calc_beta_rates(shapefacs, closed, hfb_gs=Evloss, zero_neg_str=True)

                        results_shape_factor(mgr, shapefacs, closed, write_psi=False, subdir=lab)
                        results_bare_strength(mgr, shapefacs, closed, subdir=lab)

                    # Open: rates with zeroed negstr (should already all be > 0, so this is redundant)
                    # Closed: rates from before, but if < 0 set to zero, and resum totals
                    # SKIP FOR FT_POLE: it is expected to be negative, so don't zero it
                    if not lab3:
                        results_beta_rates(mgr, shapefacs, closed, rates_df, zero_neg=True, subdir=lab)

                # For split contours, sum the contributions here
                # (FT_pole already subtracted, and rates_df will never be that of pole since its 1st in lab3_list)
                if len(split_list) > 1:
                    ratesz_df = shapefacs.zeroNegRatesBetaOut(rates_df)
                    btotz_df = btot_df.copy()
                    btotz_df.loc[btotz_df[u'B_total'] < 0, u'B_total'] = 0
                    if tot_rates_df is None:
                        tot_btot_df = btot_df.copy()
                        tot_btotz_df = btotz_df.copy()
                        tot_rates_df = rates_df.copy()
                        tot_ratesz_df = ratesz_df.copy()
                        emin = shapefacs.contour.energy_min
                        emax = shapefacs.contour.energy_max
                        all_conv = [shapefacs.all_conv]
                        total_time = shapefacs.total_time + pole_time
                    else:
                        tot_btot_df = tot_btot_df.add(btot_df)
                        tot_btotz_df = tot_btotz_df.add(btotz_df)
                        tot_rates_df = tot_rates_df.add(rates_df)
                        tot_ratesz_df = tot_ratesz_df.add(ratesz_df)
                        emin = min(emin, shapefacs.contour.energy_min)
                        emax = max(emax, shapefacs.contour.energy_max)
                        all_conv.append(shapefacs.all_conv)
                        total_time += shapefacs.total_time + pole_time

            # If we summed several split contours, write the result to directory above splits
            if tot_rates_df is not None:
                tot_rates_df = shapeFactor.recalcHalflives(tot_rates_df)
                tot_ratesz_df = shapeFactor.recalcHalflives(tot_ratesz_df)

                # Set the contour interval to the whole range for the output header
                shapefacs.contour.updateSettings({'energy_min': emin, 'energy_max': emax})
                shapefacs.sf_metadict.update({u'FAM_ctr': shapefacs.contour.name_and_int})

                results_beta_rates(mgr, shapefacs, True, tot_rates_df, zero_neg=False, subdir=ctre.label1)
                results_beta_rates(mgr, shapefacs, True, tot_ratesz_df, zero_neg=True, subdir=ctre.label1)

                tot_ct = (convString(all_conv), total_time)
                write_beta_log(mgr, shapefacs, tot_rates_df, hfb_gs, psi_log, conv_time=tot_ct, subdir=ctre.label1)

                # This code is only executed if contours are closed, so fudging is_raw is fine
                shapefacs.is_raw = True
                results_bare_strength(mgr, shapefacs, closed, btot_df=tot_btot_df, subdir=ctre.label1)
                shapefacs.is_raw = False
                results_bare_strength(mgr, shapefacs, closed, btot_df=tot_btotz_df, subdir=ctre.label1)

    except RuntimeError as rte:
        msg = ["There was an error calculating the beta decay rates."]+[str(rte)]
        pynfam_warn(msg, mgr.paths.calclabel)
    # Return ground state for dripline loop
    return hfb_gs