#!/bin/bash
#SBATCH -p pdebug
#SBATCH -J spinsurvey
#SBATCH -o out_spinsurvey_%A.out
#SBATCH -N 2
#SBATCH --ntasks-per-node=112
#SBATCH --time=00-01:00:00   #format days-hh:mm:ss

# Same conventions as slurm_ruby.sh, except this build has USE_OPENMP=1, so
# OpenMP must be explicitly disabled here -- otherwise every one of the 112
# MPI ranks' OpenMP runtime tries to grab threads for all cores on the node,
# oversubscribing it and segfaulting hfbtho_main immediately (inside libiomp5,
# __kmp_invoke_microtask). MPI-level parallelism across FAM points/blocking
# candidates is all this survey needs.
export OMP_NUM_THREADS=1

ulimit -s unlimited
srun ./run_pynfam_spinsurvey.py
