# Making useful modules available for the script
import sys
import os           # operating system module
import re           # string manipulation
import math         # math library
import collections
import itertools
import read_UNEDFdatabase01 as db

from operator import itemgetter, attrgetter

# To be set by the user:
#   - UNEDF_version .....: version of the UNEDF protocol to use. Choices are "UNEDF0", "UNEDF1", "UNEDF2", "UNEDF1-HFB"
#   - file_UNEDF ........: full name (with location) of the 'official' experimental database
#   - delta3new .........: logical flag that decides if the pairing gap is computed from odd-even mass differences or
#                          from an estimate of the HFB pairing gap (to be set to false for all choices so far)
#   - useMeas2013 .......: logical flag that decides whether to incorporate the 2013 mass measurements ("Argonne" masses) or not.  
UNEDF_version = "UNEDF1"
#file_UNEDF    = '/home/schunck/Development/sandbox/WapstraMass/data/DataSet04.dat'
file_UNEDF    = 'DataSet04.dat'
delta3new     = False
useMeas2013   = True


# Preset lists of points in the chi2: UNEDF0, UNEDF1, UNEDF2, UNEDF1-HFB

# 1) UNEDF0
Z_def_0 = [ 108, 106, 104, 102, 102, 102, 100, 100, 100, 100, 100, 100,  98,  98,  98,  98,  98, \
	     98,  98,  96,  96,  96,  94,  92,  92,  90,  72,  70,  70,  68,  68,  66,  66,  66, \
	     66,  66,  66,  66,  64,  64,  64,  64,  64,  64 ]
N_def_0 = [ 156, 154, 152, 154, 152, 150, 156, 154, 152, 150, 148, 146, 156, 154, 152, 150, 148, \
            146, 144, 150, 148, 144, 144, 144, 142, 142, 104, 108, 100, 104, 102, 102, 100,  98, \
             96,  94,  92,  90,  98,  96,  94,  92,  90,  88 ]

Z_sph_0 = [  82,  82,  82,  82,  82,  82,  82,  82,  82, 50, 50, 50, 50, 50, 50, 50, 50, 28, 28, \
	     28,  28,  28,  20,  20,  20,  20,  20,  20 ]
N_sph_0 = [ 132, 130, 128, 126, 124, 122, 120, 118, 116, 74, 72, 70, 68, 66, 64, 62, 58, 36, 34, \
             32,  30,  28,  30,  28,  26,  24,  22,  20 ]

Z_OEMn_0 = [ 100,  92,  72, 66 ]
N_OEMn_0 = [ 152, 144, 104, 98 ]

Z_OEMp_0 = [  96,  92,  68, 66 ]
N_OEMp_0 = [ 148, 142, 102, 94 ]

Z_fis_0 = [ ]
N_fis_0 = [ ]

Z_sp_0 = []
N_sp_0 = []

# 2) UNEDF1
Z_def_1 = [ 108, 106, 104, 102, 102, 102, 100, 100, 100, 100, 100, 100,  98,  98,  98,  98,  98, \
	     98,  98,  96,  96,  96,  96,  94,  94,  92,  92,  92,  90,  72,  70,  70,  68,  68, \
	     66,  66,  66,  66,  66,  66,  66,  64,  64,  64,  64,  64,  64 ]
N_def_1 = [ 156, 154, 152, 154, 152, 150, 156, 154, 152, 150, 148, 146, 156, 154, 152, 150, 148, \
            146, 144, 150, 148, 146, 144, 146, 144, 146, 144, 142, 142, 104, 108, 100, 104, 102, \
            102, 100,  98,  96,  94,  92,  90,  98,  96,  94,  92,  90,  88 ]

Z_sph_1 = [  82,  82,  82,  82,  82,  82,  82,  82,  82, 50, 50, 50, 50, 50, 50, 50, 50, 28, 28, \
	     28,  28,  28,  20,  20,  20,  20,  20,  20 ]
N_sph_1 = [ 132, 130, 128, 126, 124, 122, 120, 118, 116, 74, 72, 70, 68, 66, 64, 62, 58, 36, 34, \
             32,  30,  28,  30,  28,  26,  24,  22,  20 ]

Z_OEMn_1 = [ 100,  92,  72, 66 ]
N_OEMn_1 = [ 152, 144, 104, 98 ]

Z_OEMp_1 = [  96,  92,  68, 66 ]
N_OEMp_1 = [ 148, 142, 102, 94 ]

Z_fis_1 = [  92,  92,  94,  96 ]
N_fis_1 = [ 144, 146, 146, 146 ]

Z_2013_1 = [50, 50, 50, 52, 52, 52, 52, 54, 54, 56, 56, 56, 58, 58, 58, 62, 62]
N_2013_1 = [80, 82, 84, 82, 84, 86, 88, 84, 86, 86, 88, 90, 88, 90, 92, 96, 98]

Z_sp_1 = []
N_sp_1 = []

# 3) UNEDF2
Z_def_2 = [ 108, 106, 104, 102, 102, 102, 100, 100, 100, 100, 100, 100,  98,  98,  98,  98,  98,  \
	     98,  98,  96,  96,  96,  96,  94,  94,  92,  92,  92,  90,  72,  70,  70,  68,  68,  \
	     66,  66,  66,  66,  66,  66,  66,  64,  64,  64,  64,  64,  64 ]
N_def_2 = [ 156, 154, 152, 154, 152, 150, 156, 154, 152, 150, 148, 146, 156, 154, 152, 150, 148, \
            146, 144, 150, 148, 146, 144, 146, 144, 146, 144, 142, 142, 104, 108, 100, 104, 102, \
            102, 100,  98,  96,  94,  92,  90,  98,  96,  94,  92,  90,  88 ]

Z_sph_2 = [  82,  82,  82,  82,  82,  82,  82,  82,  82,  50, 50, 50, 50, 50, 50, 50, 50, 50, 28, \
	     28,  28,  28,  28,  20,  20,  20,  20,  20,  20 ]
N_sph_2 = [ 132, 130, 128, 126, 124, 122, 120, 118, 116,  82, 74, 72, 70, 68, 66, 64, 62, 58, 36, \
	     34,  32,  30,  28,  30,  28,  26,  24,  22,  20 ]

Z_OEMn_2 = [ 100,  90,  92,  72, 66, 50, 50 ]
N_OEMn_2 = [ 152, 142, 144, 104, 98, 74, 70 ]

Z_OEMp_2 = [  96,  92,  90, 76,  68, 66 ]
N_OEMp_2 = [ 148, 142, 142, 90, 102, 94 ]

Z_fis_2 = [  92,  92,  94,  96 ]
N_fis_2 = [ 144, 146, 146, 146 ]

Z_sp_2 = []
N_sp_2 = []


# 4) UNEDF1-HFB (identical to UNEDF1 for the time being)
Z_def_1b = [ 108, 106, 104, 102, 102, 102, 100, 100, 100, 100, 100, 100,  98,  98,  98,  98,  98, \
	      98,  98,  96,  96,  96,  96,  94,  94,  92,  92,  92,  90,  72,  70,  70,  68,  68, \
	      66,  66,  66,  66,  66,  66,  66,  64,  64,  64,  64,  64,  64 ]
N_def_1b = [ 156, 154, 152, 154, 152, 150, 156, 154, 152, 150, 148, 146, 156, 154, 152, 150, 148, \
             146, 144, 150, 148, 146, 144, 146, 144, 146, 144, 142, 142, 104, 108, 100, 104, 102, \
             102, 100,  98,  96,  94,  92,  90,  98,  96,  94,  92,  90,  88 ]

Z_sph_1b = [  82,  82,  82,  82,  82,  82,  82,  82,  82, 50, 50, 50, 50, 50, 50, 50, 50, 28, 28, \
	      28,  28,  28,  20,  20,  20,  20,  20,  20 ]
N_sph_1b = [ 132, 130, 128, 126, 124, 122, 120, 118, 116, 74, 72, 70, 68, 66, 64, 62, 58, 36, 34, \
              32,  30,  28,  30,  28,  26,  24,  22,  20 ]

Z_OEMn_1b = [ 100,  92,  72, 66 ]
N_OEMn_1b = [ 152, 144, 104, 98 ]

Z_OEMp_1b = [  96,  92,  68, 66 ]
N_OEMp_1b = [ 148, 142, 102, 94 ]

Z_fis_1b = [  92,  92,  94,  96 ]
N_fis_1b = [ 144, 146, 146, 146 ]

Z_2013_1b = [50, 50, 50, 52, 52, 52, 52, 54, 54, 56, 56, 56, 58, 58, 58, 62, 62]
N_2013_1b = [80, 82, 84, 82, 84, 86, 88, 84, 86, 86, 88, 90, 88, 90, 92, 96, 98]


Z_sp_1b = []
N_sp_1b = []


# Build dictionaries of datasets

dict_Z_def = { "UNEDF0": Z_def_0, "UNEDF1": Z_def_1, "UNEDF2": Z_def_2, "UNEDF1-HFB": Z_def_1b }
dict_N_def = { "UNEDF0": N_def_0, "UNEDF1": N_def_1, "UNEDF2": N_def_2, "UNEDF1-HFB": N_def_1b }

dict_Z_sph = { "UNEDF0": Z_sph_0, "UNEDF1": Z_sph_1, "UNEDF2": Z_sph_2, "UNEDF1-HFB": Z_sph_1b}
dict_N_sph = { "UNEDF0": N_sph_0, "UNEDF1": N_sph_1, "UNEDF2": N_sph_2, "UNEDF1-HFB": N_sph_1b}

dict_Z_OEMn = { "UNEDF0": Z_OEMn_0, "UNEDF1": Z_OEMn_1, "UNEDF2": Z_OEMn_2, "UNEDF1-HFB": Z_OEMn_1b }
dict_N_OEMn = { "UNEDF0": N_OEMn_0, "UNEDF1": N_OEMn_1, "UNEDF2": N_OEMn_2, "UNEDF1-HFB": N_OEMn_1b }

dict_Z_OEMp = { "UNEDF0": Z_OEMp_0, "UNEDF1": Z_OEMp_1, "UNEDF2": Z_OEMp_2, "UNEDF1-HFB": Z_OEMp_1b }
dict_N_OEMp = { "UNEDF0": N_OEMp_0, "UNEDF1": N_OEMp_1, "UNEDF2": N_OEMp_2, "UNEDF1-HFB": N_OEMp_1b }

dict_Z_fis = { "UNEDF0": Z_fis_0, "UNEDF1": Z_fis_1, "UNEDF2": Z_fis_2, "UNEDF1-HFB": Z_fis_1b }
dict_N_fis = { "UNEDF0": N_fis_0, "UNEDF1": N_fis_1, "UNEDF2": N_fis_2, "UNEDF1-HFB": N_fis_1b }

dict_Z_2013 = { "UNEDF1": Z_2013_1, "UNEDF1-HFB": Z_2013_1b }
dict_N_2013 = { "UNEDF1": N_2013_1, "UNEDF1-HFB": N_2013_1b }

dict_Z_sp = { "UNEDF0": Z_sp_0, "UNEDF1": Z_sp_1, "UNEDF2": Z_sp_2, "UNEDF1-HFB": Z_sp_1b }
dict_N_sp = { "UNEDF0": N_sp_0, "UNEDF1": N_sp_1, "UNEDF2": N_sp_2, "UNEDF1-HFB": N_sp_1b }

w_mass_dict    = { "UNEDF0": 2.00, "UNEDF1": 2.00, "UNEDF2": 2.00, "UNEDF1-HFB": 2.00 }
w_radius_dict  = { "UNEDF0": 0.02, "UNEDF1": 0.02, "UNEDF2": 0.02, "UNEDF1-HFB": 0.02 }
w_OEM_dict     = { "UNEDF0": 0.05, "UNEDF1": 0.05, "UNEDF2": 0.05, "UNEDF1-HFB": 0.10 }
w_fission_dict = { "UNEDF0": 0.50, "UNEDF1": 0.50, "UNEDF2": 0.50, "UNEDF1-HFB": 0.50 }
w_2013_dict    = { "UNEDF0": 2.00, "UNEDF1": 2.00, "UNEDF2": 2.00, "UNEDF1-HFB": 2.00 }
w_sp_dict      = { "UNEDF0": 1.20, "UNEDF1": 1.20, "UNEDF2": 1.20, "UNEDF1-HFB": 1.20 }


# Select dataset based on version of UNEDF fit

Z_def = dict_Z_def[UNEDF_version]
N_def = dict_N_def[UNEDF_version]

Z_sph = dict_Z_sph[UNEDF_version]
N_sph = dict_N_sph[UNEDF_version]

Z_OEMn = dict_Z_OEMn[UNEDF_version]
N_OEMn = dict_N_OEMn[UNEDF_version]

Z_OEMp = dict_Z_OEMp[UNEDF_version]
N_OEMp = dict_N_OEMp[UNEDF_version]

Z_fis = dict_Z_fis[UNEDF_version]
N_fis = dict_N_fis[UNEDF_version]

Z_2013 = dict_Z_2013[UNEDF_version]
N_2013 = dict_N_2013[UNEDF_version]

Z_sp = dict_Z_sp[UNEDF_version]
N_sp = dict_N_sp[UNEDF_version]


# Define weights in the chi2

w_mass    = w_mass_dict[UNEDF_version]
w_radius  = w_radius_dict[UNEDF_version]
w_OEM     = w_OEM_dict[UNEDF_version]
w_fission = w_fission_dict[UNEDF_version]
w_2013 = w_2013_dict[UNEDF_version]
w_sp      = w_sp_dict[UNEDF_version]


# Read experimental database

new_db = db.UNEDFdatabase(file_UNEDF)

#-------------------------------------------#
#    BUILDING THE CHI2  AND HFBTHO INPUTS   #
#-------------------------------------------#

counter = 0

# empirical conversion factor from deformation to multipole moment
#factor = (3.0*math.sqrt(2.0)*1.44)/(4*math.pi*100.0)
factor = math.sqrt(5.0/math.pi)/100.0

line_chi2 = []
line_hfbtho = []

# Deformed nuclei
dico_def = new_db.get_data("deformed")
data_defE = []
data_defb = []
for Z,N,E,b2 in zip(dico_def["Z"],dico_def["N"],dico_def["E"],dico_def["beta"]):
	key = str(Z) + "_" + str(N)
	valE = str(E)
	valb = str(b2)
	data_defE.append( (key, valE) )
	data_defb.append( (key, valb) )
dict_defE = dict(data_defE)
dict_defb = dict(data_defb)

# masses
for Z,N in zip(Z_def,N_def):
	try:
		counter=counter+1
		line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
			       + str(dict_defE[str(Z)+"_"+str(N)]) + "\t" + str(w_mass) \
			       + "\t" + "1\n")
		beta2 = float(dict_defb[str(Z)+"_"+str(N)])
		mass = float(Z)+float(N)
		key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)
		blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
		line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking ) ) )
	except KeyError:
		print "Deformed mass - Nucleus ", Z, N, " not in UNEDF database"

#  Nuclei from the 2013 mass measurement (the "Argonne" masses)
if (useMeas2013 and UNEDF_version == "UNEDF1" or UNEDF_version == "UNEDF1-HFB"):
	dico_2013 = new_db.get_data("meas2013")
	data_2013E = []
	data_2013b = []
	for Z,N,E,b2 in zip(dico_2013["Z"],dico_2013["N"],dico_2013["E"],dico_2013["beta"]):
		key = str(Z) + "_" + str(N)
		valE = str(E)
		valb = str(b2)
		data_2013E.append( (key, valE) )
		data_2013b.append( (key, valb) )
	dict_2013E = dict(data_2013E)
	dict_2013b = dict(data_2013b)

	# masses
	for Z,N in zip(Z_2013,N_2013):
		try:
			counter=counter+1
			line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
				       + str(dict_2013E[str(Z)+"_"+str(N)]) + "\t" + str(w_2013) \
				       + "\t" + "1\n")
			beta2 = float(dict_2013b[str(Z)+"_"+str(N)])
			mass = float(Z)+float(N)
			key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)
			blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
			basis_deformation_value = 0.0
			if (beta2 > 10.0):
				basis_deformation_value = 0.3
			line_hfbtho.append( ( key, (beta2, basis_deformation_value, blocking ) ) )
		except KeyError:
			print "Meas2013 mass - Nucleus ", Z, N, " not in UNEDF database"

# Spherical nuclei
dico_sph = new_db.get_data("spherical")
data_sphE = []
data_sphR = []
for Z,N,E,R in zip(dico_sph["Z"],dico_sph["N"],dico_sph["E"],dico_sph["Rp"]):
	key = str(Z) + "_" + str(N)
	valE = str(E)
	valR = str(R)
	data_sphE.append( (key, valE) )
	data_sphR.append( (key, valR) )
dict_sphE = dict(data_sphE)
dict_sphR = dict(data_sphR)

# masses
for Z,N in zip(Z_sph,N_sph):
	try:
		counter=counter+1
		line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
			       + str(dict_sphE[str(Z)+"_"+str(N)]) + "\t" + str(w_mass) \
			       + "\t" + "1\n")
		key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N)
		blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
		line_hfbtho.append( ( key, (0.0, 0.0, blocking) ) )
	except KeyError:
		print "Spherical Mass - Nucleus ", Z, N, " not in UNEDF database"

# radii
for Z,N in zip(Z_sph,N_sph):
	try:
		counter=counter+1
		line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
			       + str(dict_sphR[str(Z)+"_"+str(N)]) + "\t" + str(w_radius) \
			       + "\t" + "18\n")
	except KeyError:
		print "R.m.s. Radius - Nucleus ", Z, N, " not in UNEDF database"

# Odd-even mass differences (neutrons)
dico_d3n = new_db.get_data("delta3n")
data_d3n = []
for Z,N,E in zip(dico_d3n["Z"],dico_d3n["N"],dico_d3n["D3n"]):
	key = str(Z) + "_" + str(N)
	val = str(E)
	data_d3n.append( (key, val) )
dict_d3n = dict(data_d3n)

for Z,N in zip(Z_OEMn,N_OEMn):
	try:
		counter=counter+1
		OEMn = 0.5*(float(dict_d3n[str(Z)+"_"+str(N-1)]) + float(dict_d3n[str(Z)+"_"+str(N+1)]))
		line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
			       + str(OEMn) + "\t" + str(w_OEM) + "\t" + "12\n")
		if delta3new:
			for Np in [N-2,N,N+2]:
				mass = float(Z)+float(Np)
				beta2 = 0.3
				key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(Np) + "_n"
				blocking = [(0, 0, 0, 0, 0), (1, 0, 0, 0, 0)]
				line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking) ) )
			for Np in [N-1,N+1]:
				mass = float(Z)+float(Np)
				beta2 = 0.3
				key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(Np-1) + "_n"
				blocking = [(0, 0, 0, 0, 0), (1, 0, 0, 0, 0)]
				line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking) ) )
		else:
			if str(Z)+"_"+str(N) not in dict_defb:
				mass = float(Z)+float(N)
				beta2 = 0.3
				key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N) + "_n"
				blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
				line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking) ) )
	except KeyError:
		print "Neutron OEM - Nucleus ", Z, N, " not in UNEDF database"

# Odd-even mass differences (protons)
dico_d3p = new_db.get_data("delta3p")
data_d3p = []
for Z,N,E in zip(dico_d3p["Z"],dico_d3p["N"],dico_d3p["D3p"]):
	key = str(Z) + "_" + str(N)
	val = str(E)
	data_d3p.append( (key, val) )
dict_d3p = dict(data_d3p)

for Z,N in zip(Z_OEMp,N_OEMp):
	try:
		counter=counter+1
		OEMp = 0.5*(float(dict_d3p[str(Z-1)+"_"+str(N)]) + float(dict_d3p[str(Z+1)+"_"+str(N)]))
		line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
			       + str(OEMp) + "\t" + str(w_OEM) + "\t" + "13\n")

		if delta3new:
			for Zp in [Z-2,Z,Z+2]:
				mass = float(Zp)+float(N)
				beta2 = 0.3
				key = "{0:=03d}".format(Zp) + "_" + "{0:=03d}".format(N) + "_p"
				blocking = [(1, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
				line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking ) ) )
			for Zp in [Z-1,Z+1]:
				mass = float(Zp)+float(N)
				beta2 = 0.3
				key = "{0:=03d}".format(Zp-1) + "_" + "{0:=03d}".format(N) + "_p"
				blocking = [(1, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
				line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking ) ) )
		else:
			if str(Z)+"_"+str(N) not in dict_defb:
				mass = float(Z)+float(N)
				beta2 = 0.3
				key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N) + "_p"
				blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
				line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking ) ) )

	except KeyError:
		print "Proton OEM - Nucleus ", Z, N, " not in UNEDF database"

# Fission isomer excitation energies
if UNEDF_version in ["UNEDF1", "UNEDF2"]:

	# Fission isomer excitation energies
	dico_fis = new_db.get_data("isomers")
	data_fis = []
	for Z,N,ESD in zip(dico_fis["Z"],dico_fis["N"],dico_fis["ESD"]):
		key = str(Z) + "_" + str(N)
		val = str(ESD)
		data_fis.append( (key, val) )
	dict_fis = dict(data_fis)

	for Z,N in zip(Z_fis,N_fis):
		try:
			counter=counter+1
			line_chi2.append(str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" \
				       + str(dict_fis[str(Z)+"_"+str(N)]) + "\t" + str(w_fission) \
				       + "\t" + "1\n")
			mass = float(Z)+float(N)
			key = "{0:=03d}".format(Z) + "_" + "{0:=03d}".format(N) +"_f"
			blocking = [(0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
			# input for ground-state
			if str(Z) + "_" + str(N) not in dict_defb:
				beta2 = 0.3
				line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.3, blocking ) ) )
			# input for fission isomer
			beta2 = 0.7
			line_hfbtho.append( ( key, (factor * beta2 * (mass**(5.0/3.0)), 0.6, blocking ) ) )

		except KeyError:
			print "Fission Isomer - Nucleus ", Z, N, " not in UNEDF database"



# Write results in two files:
# -[UNEDF_VERSION]_chi2.dat contains all points in the chi2 with id, Z, N, value of simulation, weight, index of HFBTHO variable
# -[UNEDF_VERSION]_hfbtho.dat contains all inputs to HFBTHO with id, Z, N, estimated Q2, basis deformation, keys for blocking
if (useMeas2013):
	str2013 = "Argonne"
else:
	str2013 = ""
fwrite = open(UNEDF_version + str2013 + "_chi2.dat",'w')
for lignes in line_chi2:
	fwrite.write(lignes)
fwrite.close()

dict_hfbtho = dict(line_hfbtho)

fwrite = open(UNEDF_version + str2013 + "_hfbtho.dat",'w')
counter = 0
for key in sorted(dict_hfbtho):
	counter = counter+1
	broken = re.split("_",key)
	Z,N = int(broken[0]), int(broken[1])
	ligne = str(counter) + "\t" + str(Z) + "\t" + str(N) + "\t" + str(dict_hfbtho[key][0]) \
	                     + "\t" + str(dict_hfbtho[key][1]) + "\t" + str(dict_hfbtho[key][2][0]) \
	                     + "\t" + str(dict_hfbtho[key][2][1])+ "\n"
	fwrite.write(ligne)
fwrite.close()

