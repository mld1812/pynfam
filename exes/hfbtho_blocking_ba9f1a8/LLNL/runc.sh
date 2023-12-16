#!/bin/sh
# For Sierra, change OMP_NUM_THREADS=8 to OMP_NUM_THREADS=6
export OMP_NUM_THREADS=6
export KMP_STACKSIZE=1024000000 
EXECUTABLE=./anl_unedf1
srun -n 100 --ntasks-per-node=2 -c $OMP_NUM_THREADS $EXECUTABLE < /dev/null >& compute.log


