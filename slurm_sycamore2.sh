#!/bin/bash
#SBATCH -p small 
#SBATCH -J gsXXX
#SBATCH -o out_%A.out
#SBATCH -N 1 #skylake: N=2, --ntasks-per-node=40, time = 2 days.
#SBATCH --ntasks-per-node=44
#SBATCH --time=00-04:00:00   #format days-hh:mm:ss

export OMPI_MCA_coll_hcoll_enable=0
ulimit -s unlimited 
#srun python test_mpi.py 
mpirun python test_mpi.py