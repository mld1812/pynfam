#!/usr/bin/env python
import sys
sys.path.append("..")
sys.path.append("../..")
from pynfam import pynfam_residual
import numpy as np
from pounders import pounders_mpi_wrapper, stat_analysis
import pandas as pd
try:
    from mpi4py import MPI
    import_mpi = True
except ImportError:
    import_mpi = False

# code for pounders fit of GTRs; fitting variables are g0p and vpair_t0_scaled
if __name__ == '__main__':
    if import_mpi:
        comm = MPI.COMM_WORLD
        rank = comm.Get_rank()
    else:
        rank = 0

    filename = 'pynfam_input_GT.csv' # input file name
    setts = {'h0p':0.0} # fixed Landau parameters
    model = pynfam_residual(input_data=filename, override_setts_fit=setts, exes_dir="../../exes", nonConv_penalty=25)
    init_guess = (1.6, 0.0, -1.0) # initial guess of g0p, g1p and vpair_t0_scaled

    # number of observables is determined by the number of rows
    n_obs = len(pd.read_csv(filename, header=0, comment=u'#'))

    lower_bounds = (-np.inf, -np.inf, -np.inf) # lower bounds of g0p, g1p and vpair_t0_scaled
    upper_bounds = (np.inf, np.inf, 0.0) # upper bounds of g0p, g1p and vpair_t0_scaled

    # Max number of POUNDerS iterations. 
    # For a quick test, change it to a small number like 2 (also change job's time limit to a small value, like 2 hours)
    max_iter = None

    # call pounders
    fun = model.fun_g0p_g1p_vpair0
    sol, summary = pounders_mpi_wrapper(fun, init_guess, n_obs, (lower_bounds, upper_bounds), max_iter=max_iter)
    if rank == 0: print(sol, summary)

    # statistical analysis
    diff = 1e-5 # finite difference for derivatives
    stat_result = stat_analysis(fun, sol, summary['residual'], (lower_bounds, upper_bounds), diff=diff)
    if rank == 0: print(*stat_result)