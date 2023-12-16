# Making useful modules available for the script
import sys
import os           # operating system module
import re           # string manipulation
import math         # math library
import collections
import itertools

import read_UNEDFdatabase01 as db
from fit_definitions import *
from build_chi2File import *
from build_hfbthoFile import *

from operator import itemgetter, attrgetter

# To be set by the user:
#   - UNEDF_version .....: version of the UNEDF protocol to use. Choices are "UNEDF0", "UNEDF1", "UNEDF2", "UNEDF1-HFB"
#   - file_UNEDF ........: full name (with location) of the 'official' experimental database
#   - delta3new .........: logical flag that decides if the pairing gap is computed from odd-even mass differences or
#                          from an estimate of the HFB pairing gap (to be set to false for all choices so far)
#   - useMeas2013 .......: logical flag that decides whether to incorporate the 2013 mass measurements ("Argonne" masses) or not.
UNEDF_version = "UNEDF2"
file_UNEDF    = 'DataSet05.dat'
delta3new     = False
useMeas2013   = False

# Select dataset based on version of UNEDF fit
Z_def, N_def = dict_Z_def[UNEDF_version], dict_N_def[UNEDF_version]
Z_sph, N_sph = dict_Z_sph[UNEDF_version], dict_N_sph[UNEDF_version]
Z_OEMn, N_OEMn = dict_Z_OEMn[UNEDF_version], dict_N_OEMn[UNEDF_version]
Z_OEMp, N_OEMp = dict_Z_OEMp[UNEDF_version],dict_N_OEMp[UNEDF_version]
Z_fis, N_fis = dict_Z_fis[UNEDF_version], dict_N_fis[UNEDF_version]
if (useMeas2013):
	Z_2013,N_2013 = dict_Z_2013[UNEDF_version], dict_N_2013[UNEDF_version]
else:
	Z_2013,N_2013 = [], []

# Define weights for the chi2
w_mass    = w_mass_dict[UNEDF_version]
w_radius  = w_radius_dict[UNEDF_version]
w_OEM     = w_OEM_dict[UNEDF_version]
w_fission = w_fission_dict[UNEDF_version]
w_2013    = w_2013_dict[UNEDF_version]
w_sp      = w_sp_dict[UNEDF_version]

# Read experimental database
new_db = db.UNEDFdatabase(file_UNEDF)

line_hfbtho, counter = [], 1
splits_n = ( splits_exp_n, liste_splittings_n)
splits_p = ( splits_exp_p, liste_splittings_p)

#-----------------------------------------#
#   CONSTRUCT THE FILE OF HFBTHO INPUTS   #
#-----------------------------------------#
new = hfbtho_deformed(line_hfbtho, new_db, Z_def, N_def, counter, (useMeas2013, Z_2013 ,N_2013), UNEDF_version)
dico_def = new[2]
new = hfbtho_spherical(new[0], Z_sph, N_sph, new[1], UNEDF_version)
dico_sph = new[2]
new = hfbtho_OEM_n(new[0], Z_OEMn, N_OEMn, new[1], UNEDF_version, dict_defb=dico_def)
new = hfbtho_OEM_p(new[0], Z_OEMp, N_OEMp, new[1], UNEDF_version, dict_defb=dico_def)
dico_def = new[2]
if UNEDF_version in ["UNEDF1", "UNEDF2"]:
	new = hfbtho_fission(new[0], Z_fis, N_fis, new[1], UNEDF_version, dict_defb=dico_def)
	dico_def, dico_fis = new[2], new[3]
if UNEDF_version in ["UNEDF2"]:
	new = hfbtho_sp(new[0], splits_n, splits_p, new[1], UNEDF_version)
	dict_sp_n, dict_sp_p = new[2], new[3]

# Write HFBTHO input file
write_hfbtho(new[0], UNEDF_version, useMeas2013)

lines_hfbtho, dicos = extract_dicts_hfbtho(UNEDF_version, useMeas2013)

#-----------------------------------------#
#            DEFINE THE CHI2              #
#-----------------------------------------#
line_chi2, counter = [], 1
new = chi2_deformed(line_chi2, new_db, Z_def, N_def, counter, w_mass, (useMeas2013, Z_2013 ,N_2013), UNEDF_version, dicos)
new = chi2_spherical(new[0], new_db, Z_sph, N_sph, new[1], w_mass, w_radius, UNEDF_version, dicos)
new = chi2_OEM(new[0], new_db,  Z_OEMn, N_OEMn, Z_OEMp, N_OEMp, new[1], w_OEM, UNEDF_version, dicos)
if UNEDF_version in ["UNEDF1", "UNEDF2"]:
	new = chi2_fission(new[0], new_db, Z_fis, N_fis, new[1], w_fission, UNEDF_version, dicos)
if UNEDF_version in ["UNEDF2"]:
	new = chi2_sp_n(new[0], splits_n, new[1], w_sp, UNEDF_version, dict_sp_n, lines_hfbtho)
	new = chi2_sp_p(new[0], splits_p, new[1], w_sp, UNEDF_version, dict_sp_p, lines_hfbtho)

# Write chi2 file
write_chi2(new[0], UNEDF_version, useMeas2013)

