# Making useful modules available for the script
import sys
if sys.version_info < (3, 0):
	print("This parser should be run with Python3")
import os           # operating system module
import subprocess   # fine-grain control of processes
import shutil       # import shell functions
import re           # string manipulation
import math         # math library
import time
import difflib

size_node = { 'quartz': 36 }

# - computer_used ......: 'quartz'
# - bank_used ..........: 'PEM', 'fission', 'nuq', 'nucfiss', 'rmatrix'
# - queue_used .........: 'debug', 'batch'
# - disk_used ..........: 'lustre2' (LC linux)
# - number_threads .....: '9' (quartz
# - number_tasks .......: to be defined
# - number_per_node ....: '4' (quartz)
name_run        = 'hfbtho_testBench_2020-04-07/'
executable      = './hfbtho_unitTests.py'
computer_used   = 'quartz'
bank_used       = 'nucfiss'
queue_used      = 'batch'
disk_used       = 'lustre2'
number_MPI_tasks= 12
number_threads  = 9
walltime        = '00:20:00'
name_job_script = 'test.pbs'

npn = int(size_node[computer_used]/number_threads)
tcc = number_MPI_tasks*number_threads
tnc = int(math.ceil(float(tcc)/float(npn)/float(number_threads)))

number_tasks    = "{0:<d}".format(number_MPI_tasks)
number_per_node = "{0:<d}".format(npn)
total_cores     = "{0:<d}".format(tcc)
total_nodes     = "{0:<d}".format(tnc)

# Printing summary table of what the user requested
print('Name of the executable ......:', executable)
print('Computer ....................:', computer_used)
print('Bank used ...................:', bank_used)
print('Type of queue ...............:', queue_used)
print('Disk partition used .........:', disk_used)
print('MPI - Number of tasks .......:', number_tasks)
print('OpenMP - Number of threads ..:', number_threads)
print('Total number of nodes .......:', int(float(total_nodes)))
print('Total number of cores .......:', int(float(total_cores)))
print('Number of tasks per node ....:', int(float(number_per_node)))
print('Job script is ...............:', name_job_script)

# Defining the header of the batch script
header_bank = { 'PEM'    : '#SBATCH -A pbronze', \
                'fission': '#SBATCH -A fission', \
                'nucfiss': '#SBATCH -A nucfiss', \
                'nuq'    : '#SBATCH -A ndtheory', \
                'rmatrix': '#SBATCH -A nrcs' }
if computer_used in [ 'quartz' ]:
	header_queue  = { 'debug': '#SBATCH -p pdebug', 'batch': '#SBATCH -p pbatch' }
	header_output = { 'debug': '#SBATCH -o ' + computer_used + '_debug.out', \
	                  'batch': '#SBATCH -o ' + computer_used + '_batch.out' }
	header_error  = { 'debug': '#SBATCH -e ' + computer_used + '_debug.err', \
	                  'batch': '#SBATCH -e ' + computer_used + '_batch.err' }
header_walltime = { 'quartz'  : '#SBATCH -t ' + walltime }
header_nodes    = { 'quartz'  : '#SBATCH -N ' + total_nodes }

header_full = []
header_full.append('#!/bin/bash\n')
header_full.append('\n')
header_full.append(header_bank[bank_used] + '\n')
header_full.append(header_queue[queue_used] + '\n')
header_full.append(header_output[queue_used] + '\n')
header_full.append(header_error[queue_used] + '\n')
header_full.append(header_walltime[computer_used] + '\n')
header_full.append(header_nodes[computer_used] + '\n')
header_full.append('#SBATCH --export=ALL\n')
header_full.append('\n')

# Defining the command launching the job and the scratch directory where it will run
scratch_dir  = { 'quartz': '/p/' + disk_used + '/schunck1/' + name_run }
command_line = { 'quartz': 'srun -n '  + number_tasks + ' -c $OMP_NUM_THREADS $HOME/local/python-3.7.2/bin/python $EXECUTABLE ' + scratch_dir[computer_used] }
job_launch   = { 'quartz': 'sbatch ' }

# reading initialization file
with open('template_serial.txt', 'r') as fread:
	lines_template = fread.readlines()

# Adjust a few variables based on input data
chaine = 'export OMP_NUM_THREADS=\n'
position = [i for i, x in enumerate(lines_template) if x == chaine]
lines_template[position[0]] = 'export OMP_NUM_THREADS=' + str(number_threads) + '\n'

chaine = 'EXECUTABLE=\n'
position = [i for i, x in enumerate(lines_template) if x == chaine]
lines_template[position[0]] = 'EXECUTABLE=\'' + executable + '\'\n'

chaine = 'SCRATCH_DIR=\n'
position = [i for i, x in enumerate(lines_template) if x == chaine]
lines_template[position[0]] = 'SCRATCH_DIR=\'' + scratch_dir[computer_used] + '\'\n'

chaine = 'srun\n'
position = [i for i, x in enumerate(lines_template) if x == chaine]
lines_template[position[0]] = command_line[computer_used] + '\n'

# writing batch job script
fichier = 'test.pbs'
with open(fichier, 'w' ) as fwrite:
	for lines in header_full:
		fwrite.write(lines)
	for lines in lines_template:
		fwrite.write(lines)

p = subprocess.call(job_launch[computer_used] + fichier, shell=True)
print('Command line was ............:', job_launch[computer_used] + fichier)
