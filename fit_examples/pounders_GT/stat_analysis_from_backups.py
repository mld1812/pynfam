#!/usr/bin/env python
import os
import sys
sys.path.append("..")
sys.path.append("../..")
import pickle
from pynfam import pynfam_residual
import numpy as np
from pounders import stat_analysis
import pandas as pd
try:
    from mpi4py import MPI
    import_mpi = True
except ImportError:
    import_mpi = False

# For pynfam_pounders_GT, perform statistical analysis for Chi square minimization from existing solutions. 
if __name__ == '__main__':
    if import_mpi:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    else:
        rank = 0

    filename = 'pynfam_input_GT.csv' # input file name
    setts = {'h0p':0.0} # fixed Landau parameters
    model = pynfam_residual(input_data=filename, override_setts_fit=setts, exes_dir="../../exes", nonConv_penalty=25)

    lower_bounds = (-np.inf, -np.inf, -np.inf) # lower bounds of g0p, g1p and vpair_t0_scaled
    upper_bounds = (np.inf, np.inf, 0.0) # upper bounds of g0p, g1p and vpair_t0_scaled

    # statistical analysis
    if rank == 0:
        sol, res_min, _, _ = model.find_minChiSqr_in_backups(('g0p', 'g1p', 'vpair_t0_scaled'))
        if import_mpi: 
            sol, res_min = comm.bcast((sol, res_min), root=0)
        print(res_min, np.sum(res_min**2))
    else:
        sol, res_min = comm.bcast(None, root=0)
    diff = 1e-5 # finite difference for derivatives
    stat_result = stat_analysis(model.fun_g0p_g1p_vpair0, sol, res_min, (lower_bounds, upper_bounds), diff=diff)
    if rank == 0: print(*stat_result)