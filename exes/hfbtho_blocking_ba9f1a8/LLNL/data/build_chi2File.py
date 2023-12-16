# Making useful modules available for the script
import sys
import os           # operating system module
import re           # string manipulation
import math         # math library
import collections
import itertools

from fit_definitions import *
from build_hfbthoFile import *

#-------------------------------------------#
#    BUILDING THE CHI2                      #
#-------------------------------------------#

def chi2_deformed(line_chi2, new_db, Z_def, N_def, counter, w_mass, ANL2013, UNEDF_version, dicos):
	''' This method adds the data from deformed nuclei to the list of points in the chi2.
	    Inputs:
	      - line_chi2 .....: the list of points in the chi2
	      - new_db ........: an instance of the UNEDF experimental database class
	      - Z_def, N_def ..: lists of proton and neutron numbers for deformed nuclei
	      - counter .......: current position in the list of points in the chi2
	      - w_mass ........: standard deviation for masses
	      - ANL2013 .......: 3-tuple indicating if the new ANL masses must be included
	      - UNEDF_version .: version of the UNEDF functional
	      - dicos .........: tuple of dictionaries giving position of nucleus in HFBTHO inputs
	    Outputs:
	      - line_chi2 .....: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	dico_def = new_db.get_data("deformed")
	data_defE = [ (str(Z) + "_" + str(N), str(E)) for Z,N,E in zip(dico_def["Z"],dico_def["N"],dico_def["E"]) ]
	dict_defE = dict(data_defE)
	for Z,N in zip(Z_def,N_def):
		try:
			key = 'F1SPN' + "{0:>03d}".format(dicos[0]["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]) \
		            + '_F0SPN001_F0SPN001'
			line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
				       + str(dict_defE[str(Z)+"_"+str(N)]) + "\t" + str(w_mass) \
				       + "\t" + key + "\t" + "1\n")
			counter=counter+1
		except KeyError:
			print "Deformed mass - Nucleus ", Z, N, " not in UNEDF database"
	#  Nuclei from the 2013 mass measurement (the "Argonne" masses)
	if (ANL2013[0] and UNEDF_version == "UNEDF1" or UNEDF_version == "UNEDF1-HFB"):
		dico_2013 = new_db.get_data("meas2013")
		data_2013E = [ (str(Z) + "_" + str(N), str(E)) for Z,N,E in zip(dico_2013["Z"],dico_2013["N"],dico_2013["E"]) ]
		dict_2013E = dict(data_2013E)
		for Z,N in zip(ANL2013[1],ANL2013[2]):
			try:
				key = 'F1SPN' + "{0:>03d}".format(dicos[0]["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]) \
		        	    + '_F0SPN001_F0SPN001'
				line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
					       + str(dict_2013E[str(Z)+"_"+str(N)]) + "\t" + str(w_2013) \
					       + "\t" + key + "\t" + "1\n")
				counter=counter+1
			except KeyError:
				print "Meas2013 mass - Nucleus ", Z, N, " not in UNEDF database"
	return (line_chi2, counter)

def chi2_spherical(line_chi2, new_db, Z_sph, N_sph, counter, w_mass, w_radius, UNEDF_version, dicos):
	''' This method adds the data from spherical nuclei to the list of points in the chi2.
	    Inputs:
	      - line_chi2 .....: the list of points in the chi2
	      - new_db ........: an instance of the UNEDF experimental database class
	      - Z_sph, N_sph ..: lists of proton and neutron numbers for spherical nuclei
	      - counter .......: current position in the list of points in the chi2
	      - w_mass ........: standard deviation for masses
	      - w_radius ......: standard deviation for radii
	      - UNEDF_version .: version of the UNEDF functional
	      - dicos .........: tuple of dictionaries giving position of nucleus in HFBTHO inputs
	    Outputs:
	      - line_chi2 .....: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	dico_sph = new_db.get_data("spherical")
	data_sphE = [ (str(Z) + "_" + str(N), str(E)) for Z,N,E in zip(dico_sph["Z"],dico_sph["N"],dico_sph["E"])]
	dict_sphE = dict(data_sphE)
	data_sphR = [ (str(Z) + "_" + str(N), str(R)) for Z,N,R in zip(dico_sph["Z"],dico_sph["N"],dico_sph["Rp"]) if float(R) > 1.0]
	dict_sphR = dict(data_sphR)
	# masses
	for Z,N in zip(Z_sph,N_sph):
		try:
			key = 'F1SPN' + "{0:>03d}".format(dicos[1]["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]) \
		            + '_F0SPN001_F0SPN001'
			line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
				       + str(dict_sphE[str(Z)+"_"+str(N)]) + "\t" + str(w_mass) \
				       + "\t" + key + "\t" + "1\n")
			counter=counter+1
		except KeyError:
			print "Spherical Mass - Nucleus ", Z, N, " not in UNEDF database"
	# radii
	for Z,N in zip(Z_sph,N_sph):
		try:
			key = 'F1SPN' + "{0:>03d}".format(dicos[1]["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]) \
		            + '_F0SPN001_F0SPN001'
			line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
				       + str(dict_sphR[str(Z)+"_"+str(N)]) + "\t" + str(w_radius) \
				       + "\t" + key + "\t" + "18\n")
			counter=counter+1
		except KeyError:
			print "R.m.s. Radius - Nucleus ", Z, N, " not in UNEDF database"
	return (line_chi2, counter)

def chi2_OEM(line_chi2, new_db, Z_OEMn, N_OEMn, Z_OEMp, N_OEMp, counter, w_OEM, UNEDF_version, dicos):
	''' This method adds the data for odd-even mass differences to the list of points in the chi2.
	    Inputs:
	      - line_chi2 .....: the list of points in the chi2
	      - new_db ........: an instance of the UNEDF experimental database class
	      - Z_OEMn, N_OEMn : lists of proton and neutron numbers for neutron OEM
	      - Z_OEMp, N_OEMp : lists of proton and neutron numbers for proton OEM
	      - counter .......: current position in the list of points in the chi2
	      - w_OEM .........: standard deviation for all OEM
	      - UNEDF_version .: version of the UNEDF functional
	      - dicos .........: tuple of dictionaries giving position of nucleus in HFBTHO inputs
	    Outputs:
	      - line_chi2 .....: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	dico_d3n = new_db.get_data("delta3n")
	data_d3n = [ (str(Z) + "_" + str(N), str(E)) for Z,N,E in zip(dico_d3n["Z"],dico_d3n["N"],dico_d3n["D3n"])]
	dict_d3n = dict(data_d3n)
	for Z,N in zip(Z_OEMn,N_OEMn):
		found = {}
		for d in dicos:
			try:
				access = d["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]
				found = d
			except KeyError:
				pass
		if found:
			key = 'F1SPN' + "{0:>03d}".format(found["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]) \
		            + '_F0SPN001_F0SPN001'
			OEMn = 0.5*(float(dict_d3n[str(Z)+"_"+str(N-1)]) + float(dict_d3n[str(Z)+"_"+str(N+1)]))
			line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
				       + str(OEMn) + "\t" + str(w_OEM) + "\t" + key + "\t" + "12\n")
			counter=counter+1
		else:
			print "Neutron OEM - Nucleus ", Z, N, " not in UNEDF database"

	# Odd-even mass differences (protons)
	dico_d3p = new_db.get_data("delta3p")
	data_d3p = [ (str(Z) + "_" + str(N), str(E)) for Z,N,E in zip(dico_d3p["Z"],dico_d3p["N"],dico_d3p["D3p"])]
	dict_d3p = dict(data_d3p)
	for Z,N in zip(Z_OEMp,N_OEMp):
		found = {}
		for d in dicos:
			try:
				access = d["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]
				found = d
			except KeyError:
				pass
		if found:
			key = 'F1SPN' + "{0:>03d}".format(found["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]) \
		            + '_F0SPN001_F0SPN001'
			OEMp = 0.5*(float(dict_d3p[str(Z-1)+"_"+str(N)]) + float(dict_d3p[str(Z+1)+"_"+str(N)]))
			line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
				       + str(OEMp) + "\t" + str(w_OEM) + "\t" + key + "\t" + "13\n")
			counter=counter+1
		else:
			print "Proton OEM - Nucleus ", Z, N, " not in UNEDF database"
	return (line_chi2, counter)

def chi2_fission(line_chi2, new_db, Z_fis, N_fis, counter, w_fission, UNEDF_version, dicos):
	''' This method adds the data for fission isomer excitation energies to the list of points in the chi2.
	    Inputs:
	      - line_chi2 .....: the list of points in the chi2
	      - new_db ........: an instance of the UNEDF experimental database class
	      - Z_fis, N_fis ..: lists of proton and neutron numbers for fission isomers
	      - counter .......: current position in the list of points in the chi2
	      - w_fission ......: standard deviation for the fission isomers
	      - UNEDF_version .: version of the UNEDF functional
	      - dicos .........: tuple of dictionaries giving position of nucleus in HFBTHO inputs
	    Outputs:
	      - line_chi2 .....: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	if UNEDF_version in ["UNEDF1", "UNEDF2"]:
		dico_fis = new_db.get_data("isomers")
		data_fis = [ (str(Z) + "_" + str(N), str(ESD)) for Z,N,ESD in zip(dico_fis["Z"],dico_fis["N"],dico_fis["ESD"])]
		dict_fis = dict(data_fis)
		for Z,N in zip(Z_fis,N_fis):
			try:
				key = 'F1SPN' + "{0:>03d}".format(dicos[2]["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]) \
		                    +'_F1SMN' + "{0:>03d}".format(dicos[0]["{0:>03d}".format(Z) + "_" + "{0:>03d}".format(N)]) \
				    +  '_F0SPN001'
				line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
					       + str(dict_fis[str(Z)+"_"+str(N)]) + "\t" + str(w_fission) \
					       + "\t" + key + "\t" + "1\n")
				counter=counter+1
			except KeyError:
				print "Fission Isomer - Nucleus ", Z, N, " not in UNEDF database"
	return (line_chi2, counter)

def chi2_sp_n(line_chi2, splits, counter, w_sp, UNEDF_version, dict_sp, lines_hfbtho):
	''' This method adds the data for fission isomer excitation energies to the list of points in the chi2.
	    Inputs:
	      - line_chi2 .....: the list of points in the chi2
	      - counter .......: current position in the list of points in the chi2
	      - w_sp ..........: standard deviation for the fission isomers
	      - UNEDF_version .: version of the UNEDF functional
	      - dicos .........: tuple of dictionaries giving position of nucleus in HFBTHO inputs
	    Outputs:
	      - line_chi2 .....: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	df = { "pp": ("1", "1", "0"),  "hh": ("1", "1", "0"), "ph": ("1", "1", "2"), "hp": ("1", "1", "2") }
	ds = { "pp": ("P", "M", "P"),  "hh": ("M", "P", "P"), "ph": ("P", "P", "M"), "hp": ("M", "M", "P") }
	if UNEDF_version in ["UNEDF2"]:
		for nuc in splits[0]:
			Z, N, e = nuc[0], nuc[1], nuc[2]
			Z, N, block_1, block_2 = dict_sp[str(nuc)]
			c1 = format_ligne(Z, N, (0.0, 0.0, [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)], "sph"))
			c2 = format_ligne(Z, N, (0.0, 0.0, [(0, 0, 0, 0, 0), block_1[1]], "sp"))
			c3 = format_ligne(Z, N, (0.0, 0.0, [(0, 0, 0, 0, 0), block_2[1]], "sp"))
			i0 = [ i for i, x in enumerate(lines_hfbtho) if x.find(c1) > -1 ]
			i1 = [ i for i, x in enumerate(lines_hfbtho) if x.find(c2) > -1 ]
			i2 = [ i for i, x in enumerate(lines_hfbtho) if x.find(c3) > -1 ]
			typ = block_1[2] + block_2[2]
			try:
				key = 'F' + df[typ][0] + 'S' + ds[typ][0] + 'N' + "{0:>03d}".format(i1[0]+1) \
			            +'_F' + df[typ][1] + 'S' + ds[typ][1] + 'N' + "{0:>03d}".format(i2[0]+1) \
			            +'_F' + df[typ][2] + 'S' + ds[typ][2] + 'N' + "{0:>03d}".format(i0[0]+1)
				line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
					             + str(e) + "\t" + str(w_sp) + "\t" + key + "\t" + "1\n")
				counter=counter+1
			except KeyError:
				print "sp data - Nucleus ", Z, N, " not in UNEDF database"
	return (line_chi2, counter)

def chi2_sp_p(line_chi2, splits, counter, w_sp, UNEDF_version, dict_sp, lines_hfbtho):
	''' This method adds the data for fission isomer excitation energies to the list of points in the chi2.
	    Inputs:
	      - line_chi2 .....: the list of points in the chi2
	      - counter .......: current position in the list of points in the chi2
	      - w_sp ..........: standard deviation for the fission isomers
	      - UNEDF_version .: version of the UNEDF functional
	      - dicos .........: tuple of dictionaries giving position of nucleus in HFBTHO inputs
	    Outputs:
	      - line_chi2 .....: the updated list of points in the chi2
	      - counter .......: the updated position in said list
	    '''
	df = { "pp": ("1", "1", "0"),  "hh": ("1", "1", "0"), "ph": ("1", "1", "2"), "hp": ("1", "1", "2") }
	ds = { "pp": ("P", "M", "P"),  "hh": ("M", "P", "P"), "ph": ("P", "P", "M"), "hp": ("M", "M", "P") }
	if UNEDF_version in ["UNEDF2"]:
		for nuc in splits[0]:
			Z, N, e = nuc[0], nuc[1], nuc[2]
			Z, N, block_1, block_2 = dict_sp[str(nuc)]
			c1 = format_ligne(Z, N, (0.0, 0.0, [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)], "sph"))
			c2 = format_ligne(Z, N, (0.0, 0.0, [block_1[1], (0, 0, 0, 0, 0)], "sp"))
			c3 = format_ligne(Z, N, (0.0, 0.0, [block_2[1], (0, 0, 0, 0, 0)], "sp"))
			i0 = [ i for i, x in enumerate(lines_hfbtho) if x.find(c1) > -1 ]
			i1 = [ i for i, x in enumerate(lines_hfbtho) if x.find(c2) > -1 ]
			i2 = [ i for i, x in enumerate(lines_hfbtho) if x.find(c3) > -1 ]
			typ = block_1[2] + block_2[2]
			try:
				key = 'F' + df[typ][0] + 'S' + ds[typ][0] + 'N' + "{0:>03d}".format(i1[0]+1) \
			            +'_F' + df[typ][1] + 'S' + ds[typ][1] + 'N' + "{0:>03d}".format(i2[0]+1) \
			            +'_F' + df[typ][2] + 'S' + ds[typ][2] + 'N' + "{0:>03d}".format(i0[0]+1)
				line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
					       + str(e) + "\t" + str(w_sp) + "\t" + key + "\t" + "1\n")
				counter=counter+1
			except KeyError:
				print "sp data - Nucleus ", Z, N, " not in UNEDF database"
	return (line_chi2, counter)

def write_chi2(line_chi2, UNEDF_version, useMeas2013):
	''' This method writes the list of points in the chi2 in a file.
	    Inputs:
	      - line_chi2 .....: the list of points in the chi2
	      - UNEDF_version .: version of the UNEDF functional
	      - useMeas2013 ...: boolean indicating if ANL masses were included
	    '''
	if useMeas2013:
		str2013 = "Argonne"
	else:
		str2013 = ""
	fwrite = open(UNEDF_version + str2013 + "_chi2.dat",'w')
	for lignes in line_chi2:
		fwrite.write(lignes)
	fwrite.close()
	return 0



