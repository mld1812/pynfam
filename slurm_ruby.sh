#!/bin/bash
#SBATCH -p pdebug 
#SBATCH -J gsXXX
#SBATCH -o out_%A.out
#SBATCH -N 2
#SBATCH --ntasks-per-node=56
#SBATCH --time=00-01:00:00   #format days-hh:mm:ss

ulimit -s unlimited 
srun ./run_pynfam5.py

#mpirun ~/pynfam_main/exes/pnfam_main.x ~/pynfam_main/tests/150Sn_test_numcpus2/000000/fam_soln/fam_meta/P-K0/000000/P-K0.in
