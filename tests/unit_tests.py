#!/usr/bin/env python
from __future__ import unicode_literals
from pynfam.utilities.mpi_utils import do_mpi
try:
    from pynfam.utilities.mpi_utils import MPI
except:
    pass
from pynfam import pynfam_mpi_calc
from copy import deepcopy
import os
import pandas as pd
import numpy as np

rank = 0
if do_mpi: rank = MPI.COMM_WORLD.Get_rank()

"""
HFB Unit Tests:
    1. Even ground state (deformation scan)
    2. Odd ground state (blocking candidate scan)
    3. Dripline mode (2Sn dripline)
    4. Ignore non-converged hfb solution mode(s)

FAM Unit Tests:
    1. Closed contour
    2. Open contour
    3. Multiple closed contours for FT-EC (All zero, pure exp, split, pure poly, repeat)
    (4. Phase space modes)
    (5. Odd ground states)

PYNFAM Unit Tests:
    1. Load balancing
    (2. Rerun mode(s))
    (3. MPI mode)
"""

# Inputs
# -----------------------------------------------------------------------------
write_expected = False
test_subset    = None
exe_location   = "../exes"


# Define tests
# -----------------------------------------------------------------------------
unit_tests = {
'hfb1': {'directories':{'outputs': 'test_gs_even'},
         'hfb_mode': {'gs_def_scan': (-1, (-0.2,0.0,0.2))}
         },

'hfb2': {'directories':{'outputs': 'test_gs_odd'},
         'hfb_mode': {'gs_def_scan': (-2, (.2,))},
         'hfb': {'neutron_blocking': [1,0,0,0,0]}
         },

'hfb3': {'directories':{'outputs': 'test_gs_dripline'},
         'hfb_mode': {'dripline_mode':0}
         },

'hfb4': {'directories':{'outputs': 'test_gs_nonconv'},
         'hfb_mode': {'ignore_nonconv':3},
         'hfb' : {'neutron_number': (80,82,84),
                  'number_iterations': 22}
         },

'fam1': {'directories':{'outputs': 'test_ctr_closed'},
         'fam_mode': {'fam_contour': 'CIRCLE'},
         },

'fam2': {'directories':{'outputs': 'test_ctr_open'},
         'fam_mode': {'fam_contour': 'CONSTR'},
         },

'fam3': {'directories': {'outputs':         'test_ctr_ft'},
         'hfb_mode':    {'gs_def_scan':     (-1, (-0.2,0.2))},
         'fam_mode':    {'fam_contour':     'CIRCLE',
                         'beta_type':       '-'},
         'ctr':         {'nr_points':       60,
                         'energy_min':      0.0},
         'psi':         {'psi_approx':      'POLYFIT'},
         'hfb':         {'set_temperature': True,
                         'temperature':     0.5}
         },

'fam4': {'directories': {'outputs':         'test_ctr_ecft'},
         'hfb_mode':    {'gs_def_scan':     (-1, (-0.2,0.2))},
         'fam_mode':    {'fam_contour':     'CIRCLE',
                         'beta_type':       'c'},
         'ctr':         {'nr_points':       10},
         'psi':         {'psi_approx':      'POLYFIT',
                         'log(pYe)':        [5,12,15]},
         'hfb':         {'set_temperature': True,
                         'temperature':     0.5,
                         'proton_number':   36,
                         'neutron_number':  50}
         },

'py1' : {'directories':{'outputs': 'test_load_balancing'},
         'nr_parallel_calcs': None,
         'fam_mode': {'fam_contour': 'CIRCLE'},
         'ctr': {'nr_points': (10,100)}
         }
}

if test_subset:
    test_subset = np.array(test_subset, copy=False, ndmin=1)
    unit_tests = {t: unit_tests[t] for t in test_subset}

# Inputs template: Nucleus Kr116, small basis/quadrature, all
#                  modes are off/default, fam interaction is none.
# -----------------------------------------------------------------------------
template = {
    'pynfam_inputs': {
        'directories': {
            'outputs' : 'test',
            'exes'    : exe_location,
            'scratch' : './'
            },

        'nr_parallel_calcs': 1,

        'rerun_mode':0,

        'hfb_mode': {
            'gs_def_scan'    : (0,()),
            'dripline_mode'  : 0,
            'ignore_nonconv' : 0
            },

        'fam_mode': {
            'fam_contour' : None,
            'beta_type'   : '-',
            'fam_ops': "1+"
            }
       },

    'override_settings': {
        'hfb' : {'proton_number'    : 36,
                 'neutron_number'   : 80,
                 'number_of_shells' : 6,
                 'number_gauss'     : 20,
                 'number_laguerre'  : 20,
                 'number_legendre'  : 40
                 },

        'ctr' : {},

        'fam' : {'interaction_name' : 'None',
                 },

        'psi' : {},
    }
}

# -----------------------------------------------------------------------------
def update_template(d):
    """ Update template values """
    temp = deepcopy(template)
    chng_keys = d.keys()
    for k in temp.keys():    # pynfam_inputs/override_settings
        for kk in chng_keys: # keys to change
            if temp[k].get(kk) is not None:
                try:
                    temp[k][kk].update(d[kk])
                except AttributeError:
                    temp[k][kk] = d[kk]

    return temp

# -----------------------------------------------------------------------------
def compare_one_row(re, rr):
    comp = pd.DataFrame()
    skip = ['Test', 'FAM_Time', 'HFB_Time']
    success, total = 0, 0

    # Loop over results and compare
    for c in list(rr.index):
        f = False
        if c in skip:
            continue
        elif pd.isna(re[c]): # Check 1st b/c nan is also float or char
            if not pd.isna(rr[c]):
                f = True
                diff = "(NaN vs {:})".format(rr[c])
        elif isinstance(re[c], float):
            if abs(re[c] - rr[c]) > 1e-8:
                f = True
                diff = "({:.6f} vs {:.6f})".format(re[c], rr[c])
        elif isinstance(re[c], int):
            if abs(re[c] - rr[c]) != 0:
                f = True
                diff = "({:d} vs {:d})".format(re[c], rr[c])
        else:
            if re[c] != rr[c]:
                f = True
                diff = "({:} vs {:})".format(re[c], rr[c])

        if not f:
            success += 1
        else:
            comp.loc[0,c] = diff
        total += 1

    cols = list(comp.columns)
    comp.loc[0,"Test"] = rr["Test"]
    comp.loc[0, "Passed"] = str(success)+"/"+str(total)
    comp.loc[0, "Failed"] = "------>"
    comp = comp[["Test", "Passed", "Failed"]+cols]
    return comp

# -----------------------------------------------------------------------------
def check_expected(df, test_subset=None):
    expected = pd.read_csv("expected.txt", delim_whitespace=True, header=0)
    if test_subset is not None:
        tests = np.array([t in test_subset for t in expected["Test"].values])
        expected = expected[tests]
    if expected.empty: return pd.DataFrame()

    # New blocking should be string, not int
    for h in ["Z-Blk","N-Blk"]:
        expected[h] = expected[h].astype(str)

    comps = []
    for (ir, rr), (ie ,re) in zip(df.iterrows(), expected.iterrows()):
        comp = compare_one_row(re, rr)
        comps.append(comp)
    res = pd.concat(comps, axis=0, sort=False)
    res.reset_index(drop=True, inplace=True)
    return res

# -----------------------------------------------------------------------------
if __name__ == '__main__':
    logs = []
    for name, test in unit_tests.items():
        if rank == 0:
            print()
            print("{:^10}Starting Test:{:^10}{:^10}".format(" ",name," "))
            print("-"*44)
            print()
        inputs = update_template(test)
        pynfam_mpi_calc(**inputs)

        # Serialize each test
        if do_mpi: MPI.COMM_WORLD.barrier()
        if rank == 0:
            # Get outputs from master log
            outdir = inputs['pynfam_inputs']['directories']['outputs']
            logfile = os.path.join(outdir,'meta/logfile_master.dat')
            log = pd.read_csv(logfile, delim_whitespace=True, header=0)

            # Insert test name
            cols = list(log.columns)
            log['Test'] = name
            log = log[['Test']+cols]

            logs.append(log)
    if rank == 0:
        # Record results
        results = pd.concat(logs, axis=0, sort=False)
        results.reset_index(drop=True, inplace=True)
        results_str = results.to_string(header=True, index=True, col_space=3)
        with open("test_results.txt", 'w') as ex:
            ex.write(results_str+'\n')

        # Record comparison to expected
        if write_expected:
            results_str = results.to_string(header=True, index=True, col_space=3)
            with open("expected.txt", 'w') as ex:
                ex.write(results_str+'\n')
        else:
            test_results = check_expected(results, test_subset)

            test_results_str = test_results.to_string(header=True, index=True, col_space=3)
            with open("test_vs_expected.txt", 'w') as ex:
                ex.write(test_results_str+'\n')

