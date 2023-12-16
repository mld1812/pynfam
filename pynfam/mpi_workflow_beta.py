# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
# -------------- Utilities -----------------
import numpy as np
import pandas as pd
import os
from copy import deepcopy
# ----------- Relative Imports -------------
from .utilities.mpi_utils import pynfam_warn
from .outputs.ls_logger import betaLogger
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-08-28'


#-------------------------------------------------------------------------------
def check_qval(hfb_gs):
    # If the Q value is negative, dont do any calcs involving phase space
    try:
        qval = hfb_gs.soln_dict[u'HFB_Qval']
        msg  = u"Cannot compute rates for Q < 0! Setting phase space to zero."
    except KeyError:
        qval = -999.9
        msg  = u"No Q value was found in the HFB solution. Phase space will be set to zero."
    if qval < 0 and not (hfb_gs.beta == u'c' and hfb_gs.ft_active):
        pynfam_warn(msg)

#-------------------------------------------------------------------------------
def results_shape_factor(mgr, shapefacs, ctr_closed, write_psi=True, subdir=u''):

    # Set the filenames and dest
    if shapefacs.is_raw or ctr_closed:
        sf_dest = os.path.join(mgr.paths.beta, subdir)
        prefix = u''
    else:
        sf_dest = os.path.join(mgr.paths.beta_m, subdir)
        prefix = u'zeroed_'

    # Files that always go to meta
    meta_dest = os.path.join(mgr.paths.beta_m, subdir)

    # Store the pswsf (real part as meta data)
    title = u"# Nuclear Beta Decay Phase Space Weighted Shape Factors (Imaginary)"
    fname = prefix+u"shapefactor_im.out"
    shapefacs.writeOutput(shapefacs.sf_df_imag, title, fname, sf_dest)

    title = u"# Nuclear Beta Decay Phase Space Weighted Shape Factors (Real)"
    fname = prefix+u"shapefactor_re.out"
    shapefacs.writeOutput(shapefacs.sf_df_real, title, fname, meta_dest)

    # Store the phase space (meta data)
    if write_psi:
        title = u"# Nuclear Beta Decay Phase Space (Imaginary)"
        fname = u"phasespace_im.out"
        shapefacs.writeOutput(shapefacs.ps_df_imag, title, fname, meta_dest)

        title = u"# Nuclear Beta Decay Phase Space (Real)"
        fname = u"phasespace_re.out"
        shapefacs.writeOutput(shapefacs.ps_df_real, title, fname, meta_dest)

        title = u"# Nuclear Beta Decay Integrated Phase Space"
        fname = u"phasespace_int.out"
        ps_int = shapefacs.integratePhaseSpace()
        shapefacs.writeOutput(ps_int, title, fname, meta_dest)

#-------------------------------------------------------------------------------
def calc_beta_rates(shapefacs, ctr_closed, hfb_gs=None, zero_neg_str=False):

    rates_df = shapefacs.calcBetaRates()

    if hfb_gs is not None:
        # We recalculate the shapefactor here, so don't zero if closed contour
        sf_evloss = deepcopy(shapefacs)
        zero_evloss = zero_neg_str
        if ctr_closed and zero_neg_str:
            zero_evloss = False
        sf_evloss.calcShapeFactor(hfb_gs, Evloss=True, zero_neg=zero_evloss)

        # Now get the rates
        evloss_df = sf_evloss.calcBetaRates()

        # Combine Rates and Evloss into one dataframe
        newnames = {h:"Evloss-"+h for h in evloss_df.columns}
        evloss_df.rename(columns=newnames, inplace=True)
        rates_df = pd.concat([rates_df, evloss_df], axis=1)

    return rates_df

#-------------------------------------------------------------------------------
def results_beta_rates(mgr, shapefacs, ctr_closed, rates_df, zero_neg=False, subdir=u''):

    # Set the filenames and dest
    if not zero_neg:
        sf_dest = os.path.join(mgr.paths.beta, subdir)
        prefix = u''
    else:
        sf_dest = os.path.join(mgr.paths.beta_m, subdir)
        prefix = u'zeroed_'

    # Store the rates
    title = u"# Nuclear Beta Decay Rates and Half-Lives"
    fname = prefix+u"beta.out"
    if zero_neg:
        rates_df = shapefacs.zeroNegRatesBetaOut(rates_df)
    shapefacs.writeOutput(rates_df, title, fname, sf_dest)

    # Store the cumulative rates
    if not ctr_closed:
        title = u"# Nuclear Beta Decay Cumulative Rates"
        fname = prefix+u"cumulative_rate.out"
        crate = shapefacs.calcCumulativeStr(rate=True)
        shapefacs.writeOutput(crate, title, fname, sf_dest)

    return rates_df

#-------------------------------------------------------------------------------
def results_bare_strength(mgr, shapefacs, ctr_closed, btot_df=None, subdir=u''):
    """ Calculate total/cumulative GT strength without phase space """

    # Set the filenames and dest
    if shapefacs.is_raw:
        sf_dest = os.path.join(mgr.paths.beta, subdir)
        prefix = u''
        zero_neg = False
    else:
        sf_dest = os.path.join(mgr.paths.beta_m, subdir)
        prefix = u'zeroed_'
        zero_neg = True
    genops = shapefacs.genops

    # Bare Gamow-Teller strength (GTK0 + 2GTK1). Writes BGT vs E and sum(BGT) vs E
    if not ctr_closed and (u'GT_K0' in genops and u'GT_K1' in genops):
        gt_str = shapefacs.calcTotalGT(zero_neg=zero_neg)
        # Update E 1st peak to be for GT, not PSWSF, for file header
        if shapefacs.is_raw:
            ind = shapefacs.findFirstPeak(gt_str[u'Re(EQRPA)'].values, gt_str[u'Total-GT'].values)[0]
            shapefacs.sf_metadict[u'E_1stPeak'] = gt_str[u'Re(EQRPA)'].values[ind]

        title = u"# Nuclear Beta Decay Gamow-Teller Strength"
        fname = prefix+u"gamow_teller.out"
        shapefacs.writeOutput(gt_str, title, fname, sf_dest)

    # Integrate bare strengths (gamow-teller info for closed is in here)
    if btot_df is None:
        btot_df = shapefacs.calcTotalStr(zero_neg=zero_neg)
    title = u"# Nuclear Beta Decay Total Strengths (PRC 90 024308 (2014) A2)"
    fname = prefix+u"total_str.out"
    shapefacs.writeOutput(btot_df, title, fname, sf_dest)


#-------------------------------------------------------------------------------
def write_beta_log(mgr, shapefacs, rates_df, hfb_gs, psi_log={}, subdir=u'', conv_time=None):
    betalog = betaLogger(u"logfile_mini.dat")
    beta_soln_dict = hfb_gs.soln_dict.copy()

    # Ignore the 0/0 warning for %FF
    with np.errstate(divide=u'ignore', invalid=u'ignore'):
        beta_soln_dict.update({
            u'Total-HL' : u'{:.6e}'.format(rates_df.loc[u'Total', u'Half-Life(s)']),
            u'%FF'      : (rates_df.loc[u'Total-Forbidden',u'Rate(s^-1)']/rates_df.loc[u'Total', u'Rate(s^-1)'])*100.0,
            u'FAM_Qval' : shapefacs.sf_metadict[u'FAM_Qval']
            })

    # Add any extra params we want
    beta_soln_dict.update(psi_log)

    # Get conv and time for all operators
    if conv_time is None:
        conv = shapefacs.all_conv
        time = shapefacs.total_time
    else:
        conv, time = conv_time
    beta_soln_dict[u'FAM_Conv'] = conv
    beta_soln_dict[u'FAM_Time'] = time

    # Update the hfb label to be that of the beta_soln
    beta_soln_dict[u'Label'] = mgr.paths.rp(os.path.join(mgr.paths.beta, subdir))

    betalog.quickWrite(beta_soln_dict, dest=os.path.join(mgr.paths.beta, subdir))
