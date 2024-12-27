#!/bin/bash
#SBATCH -p small
#SBATCH -J gsXXX
#SBATCH -o out_%A.out
#SBATCH -N 1 #skylake: N=2, --ntasks-per-node=40, time = 2 days.
#SBATCH --ntasks-per-node=26
##SBATCH --cpus-per-task=1
#SBATCH --time=00-04:00:00   #format days-hh:mm:ss
export OMPI_MCA_coll_hcoll_enable=0

ulimit -s unlimited 
mpirun ./run_pynfam5.py

#mpirun ~/pynfam_main/exes/pnfam_main.x ~/pynfam_main/tests/150Sn_test_numcpus2/000000/fam_soln/fam_meta/P-K0/000000/P-K0.in
