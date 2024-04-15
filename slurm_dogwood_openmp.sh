#!/bin/bash
#SBATCH -p debug_queue 
#SBATCH -J gsXXX
#SBATCH -o out_%A.out
#SBATCH --time=00-04:00:00   #format days-hh:mm:ss
#SBATCH -N 2     #previous tests (601, 605): N=1, ntasks per node = 1. old test with this one (605): n=2, tasks/node=2, cpus/task = 22, num threads = 44.
#SBATCH --ntasks-per-node=44
#SBATCH --cpus-per-task=1
##SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread

#export OMP_NUM_THREADS=4
#export OMP_STACKSIZE=256m
#export OMP_DYNAMIC=false
#export OMP_PROC_BIND=true
#export OMP_WAIT_POLICY=active

ulimit -s unlimited
mpirun ./run_pynfam2.py
