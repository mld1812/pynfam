# Making useful modules available for the script
import sys
import os           # operating system module
import re           # string manipulation
import math         # math library
import collections
import itertools
import read_UNEDFdatabase01 as db

from operator import itemgetter, attrgetter

#-------------------------------------------#
#    BUILDING THE CHI2  AND HFBTHO INPUTS   #
#-------------------------------------------#

def hfbtho_deformed(line_hfbtho, new_db, Z_def, N_def, counter, ANL2013, UNEDF_version):
	''' This method adds deformed nuclei to the list of HFBTHO inputs
	    Inputs:
	      - line_hfbtho ...: the list of points to run by HFBTHO
	      - new_db ........: an instance of the UNEDF experimental database class (needed to extract
	                         the deformation of the current nucleus)
	      - Z_def, N_def ..: lists of proton and neutron numbers for deformed nuclei
	      - counter .......: current position in the HFBTHO input list
	      - ANL2013 .......: 3-tuple indicating if the new ANL masses must be included
	      - UNEDF_version .: version of the UNEDF functional
	    Outputs:
	      - line_hfbtho ...: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	# empirical conversion factor from deformation to multipole moment
	#factor = (3.0*math.sqrt(2.0)*1.44)/(4*math.pi*100.0)
	factor = math.sqrt(5.0/math.pi)/100.0
	# Deformed nuclei
	dico_def = new_db.get_data("deformed")
	data_defb = [ (str(Z) + "_" + str(N), str(b2)) for Z,N,b2 in zip(dico_def["Z"],dico_def["N"],dico_def["beta"]) ]
	dict_defb = dict(data_defb)
	dico_def = {}
	for Z,N in zip(Z_def,N_def):
		try:
			beta2 = float(dict_defb[str(Z)+"_"+str(N)])
			mass = float(Z)+float(N)
			key = "{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)
			blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
			line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking, "def" ) ) )
			dico_def[key] = counter
			counter=counter+1
		except KeyError:
			print "Deformed mass - Nucleus ", Z, N, " not in UNEDF database"

	#  Nuclei from the 2013 mass measurement (the "Argonne" masses)
	if (ANL2013[0] and UNEDF_version == "UNEDF1" or UNEDF_version == "UNEDF1-HFB"):
		dico_2013 = new_db.get_data("meas2013")
		data_2013b = [ (str(Z) + "_" + str(N), str(b2)) for Z,N,b2 in zip(dico_2013["Z"],dico_2013["N"],dico_2013["beta"])]
		dict_2013b = dict(data_2013b)
		for Z,N in zip(ANL2013[1],ANL2013[2]):
			try:
				beta2 = float(dict_2013b[str(Z)+"_"+str(N)])
				mass = float(Z)+float(N)
				key = "{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)
				blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
				line_hfbtho.append( ( key, (beta2, 0.3, blocking, "def" ) ) )
				dico_def[key] = counter
				counter=counter+1
			except KeyError:
				print "Meas2013 mass - Nucleus ", Z, N, " not in UNEDF database"
	return (line_hfbtho, counter, dico_def)


def hfbtho_spherical(line_hfbtho, Z_sph, N_sph, counter, UNEDF_version):
	''' This method adds spherical nuclei to the list of HFBTHO inputs
	    Inputs:
	      - line_hfbtho ...: the list of points to run by HFBTHO
	      - Z_sph, N_sph ..: lists of proton and neutron numbers for spherical nuclei
	      - counter .......: current position in the HFBTHO input list
	      - UNEDF_version .: version of the UNEDF functional
	    Outputs:
	      - line_hfbtho ...: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	dico_sph = {}
	for Z,N in zip(Z_sph,N_sph):
		try:
			key = "{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)
			blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
			line_hfbtho.append( ( key, (0.0, 0.0, blocking, "sph") ) )
			dico_sph[key] = counter
			counter=counter+1
		except KeyError:
			print "Spherical Mass - Nucleus ", Z, N, " not in UNEDF database"
	return (line_hfbtho, counter, dico_sph)


def hfbtho_OEM_n(line_hfbtho, Z_OEMn, N_OEMn, counter, UNEDF_version, dict_defb={}, delta3new=False):
	''' This method adds neutron OEM to the list of HFBTHO inputs
	    Inputs:
	      - line_hfbtho ...: the list of points to run by HFBTHO
	      - Z_OEMn, N_OEMn : lists of proton and neutron numbers for neutron OEM
	      - counter .......: current position in the HFBTHO input list
	      - UNEDF_version .: version of the UNEDF functional
	      - dict_defb .....: dictionary indicating if these nuclei have already been added to the list
	      - delta3new .....: boolean indicating if theoretical OEM should be computed from 3 point formulas
	    Outputs:
	      - line_hfbtho ...: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	factor = math.sqrt(5.0/math.pi)/100.0
	for Z,N in zip(Z_OEMn,N_OEMn):
		try:
			if delta3new:
				for Np in [N-2,N,N+2]:
					mass = float(Z)+float(Np)
					key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(Np) + "_n"
					blocking = [(0, 0, 0, 0, 0), (1, 0, 0, 0, 0)]
					line_hfbtho.append( ( key, (factor * 0.3 * (mass**(5.0/3.0)), 0.3, blocking, "OEMn") ) )
					dict_defb["{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)] = counter
					counter=counter+1
				for Np in [N-1,N+1]:
					mass = float(Z)+float(Np)
					key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(Np-1) + "_n"
					blocking = [(0, 0, 0, 0, 0), (1, 0, 0, 0, 0)]
					line_hfbtho.append( ( key, (factor * 0.3 * (mass**(5.0/3.0)), 0.3, blocking, "OEMn") ) )
					dict_defb["{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)] = counter
					counter=counter+1
			else:
				if "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N) not in dict_defb:
					mass = float(Z)+float(N)
					key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N) + "_n"
					blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
					line_hfbtho.append( ( key, (factor * 0.3 * (mass**(5.0/3.0)), 0.3, blocking, "OEMn") ) )
					dict_defb["{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)] = counter
					counter=counter+1
		except KeyError:
			print "Neutron OEM - Nucleus ", Z, N, " not in UNEDF database"
	return (line_hfbtho, counter, dict_defb)


def hfbtho_OEM_p(line_hfbtho, Z_OEMp, N_OEMp, counter, UNEDF_version, dict_defb={}, delta3new=False):
	''' This method adds proton OEM to the list of HFBTHO inputs
	    Inputs:
	      - line_hfbtho ...: the list of points to run by HFBTHO
	      - Z_OEMp, N_OEMp : lists of proton and neutron numbers for proton OEM
	      - counter .......: current position in the HFBTHO input list
	      - UNEDF_version .: version of the UNEDF functional
	      - dict_defb .....: dictionary indicating if these nuclei have already been added to the list
	      - delta3new .....: boolean indicating if theoretical OEM should be computed from 3 point formulas
	    Outputs:
	      - line_hfbtho ...: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	factor = math.sqrt(5.0/math.pi)/100.0
	for Z,N in zip(Z_OEMp,N_OEMp):
		try:
			if delta3new:
				for Zp in [Z-2,Z,Z+2]:
					mass = float(Zp)+float(N)
					key = "{0:=03d}".format(Zp) + "_" + "{0:=03d}".format(N) + "_p"
					blocking = [(1, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
					line_hfbtho.append( ( key, (factor * 0.3 * (mass**(5.0/3.0)), 0.3, blocking, "OEMp" ) ) )
					dict_defb["{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)] = counter
					counter=counter+1
				for Zp in [Z-1,Z+1]:
					mass = float(Zp)+float(N)
					key = "{0:=03d}".format(Zp-1) + "_" + "{0:=03d}".format(N) + "_p"
					blocking = [(1, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
					line_hfbtho.append( ( key, (factor * 0.3 * (mass**(5.0/3.0)), 0.3, blocking, "OEMp" ) ) )
					dict_defb["{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)] = counter
					counter=counter+1
			else:
				if "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N) not in dict_defb:
					mass = float(Z)+float(N)
					key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N) + "_p"
					blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
					line_hfbtho.append( ( key, (factor * 0.3 * (mass**(5.0/3.0)), 0.3, blocking, "OEMp" ) ) )
					dict_defb["{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)] = counter
					counter=counter+1
		except KeyError:
			print "Proton OEM - Nucleus ", Z, N, " not in UNEDF database"
	return (line_hfbtho, counter, dict_defb)


def hfbtho_fission(line_hfbtho, Z_fis, N_fis, counter, UNEDF_version, dict_defb={}):
	''' This method adds fission isomers to the list of HFBTHO inputs
	    Inputs:
	      - line_hfbtho ...: the list of points to run by HFBTHO
	      - Z_fis, N_fis ..: lists of proton and neutron numbers for fission isomers
	      - counter .......: current position in the HFBTHO input list
	      - UNEDF_version .: version of the UNEDF functional
	      - dict_defb .....: dictionary indicating if these nuclei have already been added to the list
	    Outputs:
	      - line_hfbtho ...: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	factor = math.sqrt(5.0/math.pi)/100.0
	dict_fiss = {}
	if UNEDF_version in ["UNEDF1", "UNEDF2"]:
		for Z,N in zip(Z_fis,N_fis):
			try:
				mass = float(Z)+float(N)
				key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N) +"_f"
				blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
				# input for ground-state
				if str(Z) + "_" + str(N) not in dict_defb:
					line_hfbtho.append( ( "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N), (factor * 0.3 * (mass**(5.0/3.0)), 0.3, blocking, "def" ) ) )
					dict_defb["{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)] = counter
					counter=counter+1
				# input for fission isomer
				line_hfbtho.append( ( key, (factor * 0.7 * (mass**(5.0/3.0)), 0.6, blocking, "fiss" ) ) )
				dict_fiss["{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)] = counter
				counter=counter+1
			except KeyError:
				print "Fission Isomer - Nucleus ", Z, N, " not in UNEDF database"
	return (line_hfbtho, counter, dict_defb, dict_fiss)

def hfbtho_sp(line_hfbtho, splits_n, splits_p, counter, UNEDF_version):
	''' This method adds single-particle data to the list of HFBTHO inputs
	    Inputs:
	      - line_hfbtho ...: the list of points to run by HFBTHO
	      - splits ........: 2-tuple containing all information about s.p. splittings
	      - counter .......: current position in the HFBTHO input list
	      - UNEDF_version .: version of the UNEDF functional
	    Outputs:
	      - line_hfbtho ...: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	if UNEDF_version in ["UNEDF2"]:
		# Neutrons
		n, list_blocks, dict_sp_n = 0, [], {}
		for nuc, conf in zip(splits_n[0], splits_n[1]):
			Z, N, e = nuc[0], nuc[1], nuc[2]
			block_1, block_2 = conf[0], conf[1]
			key = "{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)
			if block_1 not in list_blocks:
				n = n+1
				line_hfbtho.append( ( key + '_' + str(n), (0.0, 0.0, [(0, 0, 0, 0, 0), block_1[1]], "sp") ) )
				list_blocks.append(block_1)
				counter=counter+1
			if block_2 not in list_blocks:
				n = n+1
				line_hfbtho.append( ( key + '_' + str(n), (0.0, 0.0, [(0, 0, 0, 0, 0), block_2[1]], "sp") ) )
				list_blocks.append(block_2)
				counter=counter+1
			dict_sp_n[str(nuc)] = (Z, N, block_1, block_2)
		# Protons
		list_blocks, dict_sp_p = [], {}
		for nuc, conf in zip(splits_p[0], splits_p[1]):
			Z, N, e = nuc[0], nuc[1], nuc[2]
			block_1, block_2 = conf[0], conf[1]
			key = "{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)
			if block_1 not in list_blocks:
				n = n+1
				line_hfbtho.append( ( key + '_' + str(n), (0.0, 0.0, [block_1[1], (0, 0, 0, 0, 0)], "sp") ) )
				list_blocks.append(block_1)
				counter=counter+1
			if block_2 not in list_blocks:
				n = n+1
				line_hfbtho.append( ( key + '_' + str(n), (0.0, 0.0, [block_2[1], (0, 0, 0, 0, 0)], "sp") ) )
				list_blocks.append(block_2)
				counter=counter+1
			dict_sp_p[str(nuc)] = (Z, N, block_1, block_2)
	return (line_hfbtho, counter, dict_sp_n, dict_sp_p)

def write_hfbtho(line_hfbtho, UNEDF_version, useMeas2013):
	''' This method writes the list of HFBTHO input points
	    Inputs:
	      - line_hfbtho ...: the list of HFBTHO input points
	      - UNEDF_version .: version of the UNEDF functional
	      - useMeas2013 ...: boolean indicating if ANL masses were included
	    '''
	if (useMeas2013):
		str2013 = "Argonne"
	else:
		str2013 = ""
	dict_hfbtho = dict(line_hfbtho)
	fwrite = open(UNEDF_version + str2013 + "_hfbtho.dat",'w')
	counter = 0
	for key in sorted(dict_hfbtho):
		counter = counter+1
		broken = re.split("_",key)
		Z,N = int(broken[0]), int(broken[1])
		ligne = str(counter) + "\t" + format_ligne(Z,N,dict_hfbtho[key])  + "\n"
		fwrite.write(ligne)
	fwrite.close()

def format_ligne(Z,N,t):
	''' Return a string from a tuple with the format of a HFBTHO input'''
	return str(Z) + "\t" + str(N) + "\t" + str(t[0]) + "\t" + str(t[1]) + "\t" + str(t[2][0]) + "\t" + str(t[2][1])+ "\t" + str(t[3])

def extract_dicts_hfbtho(UNEDF_version, useMeas2013):
	''' This method writes the list of HFBTHO input points
	    Inputs:
	      - UNEDF_version .: version of the UNEDF functional
	      - useMeas2013 ...: boolean indicating if ANL masses were included
	    '''
	if (useMeas2013):
		str2013 = "Argonne"
	else:
		str2013 = ""
	fread = open(UNEDF_version + str2013 + "_hfbtho.dat",'r')
	all_lines = fread.readlines()
	fread.close()

	counter = 0
	dico_sph, dico_def, dico_fis = {}, {}, {}
	for line in all_lines:
		counter = counter+1
		broken = re.split("\t",re.split("\n",line)[0])
		num, Z,N, typ = int(broken[0]), int(broken[1]), int(broken[2]), broken[7]
		if typ == 'sph':
			dico_sph["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)] = num
		if typ == 'fiss':
			dico_fis["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)] = num
		if typ in ['def', 'OEMn', 'OEMp']:
			dico_def["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)] = num
		dicos = (dico_def, dico_sph, dico_fis)
	return (all_lines, dicos)

