# Making useful modules available for the script
import sys
if sys.version_info < (3, 0):
	print("This parser should be run with Python3")
import os           # operating system module
import subprocess   # fine-grain control of processes
import shutil       # import shell functions
import re           # string manipulation
import itertools
import time
import glob
from collections import defaultdict,deque
try:
	from mpi4py import MPI
	do_mpi = True
except:
	do_mpi = False

HOME = os.path.expanduser("~")

directory_inputs  = HOME + '/git/hfbtho/src/unit_tests/'
directory_execs   = HOME + '/git/hfbtho/src/hfbtho/'
directory_outputs = directory_inputs
run_executable    = False
analyze_run       = True

liste_fichiers = os.listdir(directory_inputs)
liste_reelle   = list(filter(lambda x: x.find("hfbtho_NAMELIST") > -1, liste_fichiers))
liste_reelle.sort()
# List of HFBTHO executables to test
liste_exec = [ 'hfbtho_main' ]
# Full list of tasks
all_tasks = itertools.product(liste_exec, [ fichier_courant for fichier_courant in liste_reelle ])

# ------------------------------ #
#  EXECUTION OF A LIST OF TASKS  #
# ------------------------------ #
list_tasks, list_procs = [], []
counter = 0
if run_executable:
	for task in all_tasks:
		counter=counter+1
		# Setup condition for parallel tasks
		condition = False
		if do_mpi:
			comm = MPI.COMM_WORLD
			rank = comm.Get_rank()
			size = comm.Get_size()
			if rank % size == (counter-1) % size:
				condition = True
		else:
			condition = True
		# Prepare current input file and run it
		if condition:
			directory = './run_' + "{0:>06d}".format(counter) + '/'
			new_exec, new_data, new_output = task[0], task[1], re.split("\.d",task[1])[0] + '.out'
			if os.path.exists(directory) == False:
				os.mkdir(directory)
			shutil.copy(directory_execs + new_exec, directory)
			shutil.copy(directory_inputs + new_data, directory + 'hfbtho_NAMELIST.dat')
			os.chdir(directory)
			# Run HFBTHO on current file
			hfbtho_argument = new_exec + " > " + new_output
			start_time = time.time()
			completed = subprocess.run(hfbtho_argument, shell=True)
			end_time = time.time()
			if completed.returncode == 0:
				print("Process ", counter, " for input file ", new_data, " and executable ", new_exec, " has completed in ", end_time - start_time)
			else:
				print("Process ", counter, " for input file ", new_data, " and executable ", new_exec, " has failed with error code ", completed.returncode)
			os.chdir("./..")

# ------------------------------ #
#  ANALYSING THE UNIT TESTS      #
# ------------------------------ #

# Function used to strip a line from its blank spaces
def not_empty(chaine): return chaine != ''

dico_obs = { "Nshell": ("Number of HO shells ........:", 5) , \
             "b0"    : ("HO length b0 (fm) ..........:", 5) , \
             "beta2" : ("Basis deformation ..........:", 4) , \
             "EDF"   : ("Energy functional ..........:", 3) , \
             "N"     : ("Requested part.numbs.", 2) , \
             "Z"     : ("Requested part.numbs.", 3) , \
             "deltan": ("delta(n,p), pwi .....", 3) , \
             "deltap": ("delta(n,p), pwi .....", 4) , \
             "Ecut"  : ("delta(n,p), pwi .....", 5) , \
             "T"     : ("Temperature", 3) , \
             "Sn"    : ("Entropy", 2) , \
             "Sp"    : ("Entropy", 3) , \
             "q2"    : ("quadrupole moment[b]", 4) , \
             "EHFB"  : ("tEnergy: ehfb (qp)...", 3) }

liste_obs = [ 'Z', 'N', 'q2', 'EHFB', 'T', 'Sn', 'Sp']
liste_format = [ "{0:>6.2f}", "{0:>6.2f}", "{0:>7.3f}", "{0:>12.6f}", "{0:>4.2f}", "{0:>6.2f}","{0:>6.2f}" ]
dico = defaultdict(list)

# ------------------------------ #
#  ANALYSIS OF THE RESULTS       #
# ------------------------------ #
if analyze_run:
	os.chdir(directory_outputs)
	lignes = []
	for rep in sorted(glob.glob('run_00*')):
		os.chdir(rep)
		fichier = sorted(glob.glob('*.out'))[0]
		# Open current file
		fread = open(fichier,'r')
		allLines = fread.readlines()
		fread.close()
		# Extract a number of useful quantities
		for key in liste_obs:
			key_obs = dico_obs[key][0]
			pos_obs = dico_obs[key][1]
			pos = [ i for i, x in enumerate(allLines) if x.find(key_obs) > -1 ]
			size = len(pos)
			deque((list.pop(allLines, i) for i in sorted(pos[0:size-1], reverse=True)), maxlen=0)
			pos = [ i for i, x in enumerate(allLines) if x.find(key_obs) > -1 ]
			if size > 0:
				currentLine  = re.split("\n",allLines[pos[0]])
				brokenLine   = re.split(' ',currentLine[0])
				strippedLine = list(filter(not_empty, brokenLine))
				dico[key].append(strippedLine[pos_obs])
		os.chdir("./..")
	# Format the outputs
	lignes = [ '#  Z      N      q2        EHFB     T      Sn     Sp  \n']
	all_listes = [ dico[k] for k in liste_obs ]
	all_obs = zip(*all_listes)
	for line in all_obs:
		line_str = [ f.format(float(field)) for f,field in zip(liste_format,line) ]
		#line_str = [ field for field in line ]
		lignes.append( ' '.join(line_str) + '\n')
	# Write new file on disk
	fwrite = open('results.txt','w')
	for l in lignes:
		fwrite.write(l)
	fwrite.close()

