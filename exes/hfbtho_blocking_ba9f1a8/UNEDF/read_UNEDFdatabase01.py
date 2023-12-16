#---------------------------------------------------------------#
#                                                               #
#              DESCRIPTION OF THE PROGRAM                       #
#                                                               #
#  Language      : Python >= 2.6                                #
#  Author        : N. Schunck                                   #
#  Last Modified : February 11, 2013                            #
#  Version       : 1.0                                          #
#  Description   :
#                                                               #
#---------------------------------------------------------------#

# Making useful modules available for the script
import sys
import os           # operating system module
import re           # string manipulation
import math         # math library
import collections
import itertools

from operator import itemgetter, attrgetter

from numpy import * # Array functions

# Function which strips a character string of its blank spaces
def not_empty(chaine): return chaine != ''

class UNEDFdatabase:

	"""The official UNEDF experimental database.

	   Attributes:
	       fichier 		name of the file that contains the dataset (with full path)

	   Methods:
	       read_file	reads the file with the dataset
	       get_data		wrapper function that selects and returns the kind of data needed
	       get_spherical	returns tuple of lists for data in spherical nuclei
               get_deformed	returns tuple of lists for data in deformed nuclei
	       get_delta3n	returns tuple of lists for neutron odd-even mass differences
	       get_delta3p	returns tuple of lists for proton odd-even mass differences
	       get_meas2013     returns tuple of lists for the 2013 mass measurements (the "Argonne" masses)

	"""
	def __init__(self, fichier):
		# inputs
		self.fichier   = fichier
		self.separator = ' '
		# outputs
		self.liste_data = []
		self.dict_data  = {}


	# Open the file, read its content and put the results into a big list
	def read_file(self):

		fichier = self.fichier
		fread = open( fichier )
		self.liste_data = fread.readlines()
		fread.close


	# Based on keyword, return data
	def get_data(self, keyword):
		self.read_file()
		internal_dict = { "spherical": self.get_spherical, "deformed": self.get_deformed, \
		                  "delta3n": self.get_delta3n, "delta3p": self.get_delta3p, \
		                  "isomers": self.get_isomers, \
				  "meas2013": self.get_meas2013 }
		return internal_dict[keyword]()


	# Obtain data for spherical nuclei
	def get_spherical(self):
		chaine = "SPHERICA !   Z   N     B [MeV]     dB[MeV]    R0 [fm]   sigma [fm]  Rch [fm]   Rp [fm]"
		position = [ i for i, x in enumerate(self.liste_data) if x.find(chaine) > -1 ]
		debut = int(position[0])+2
		N = int(self.liste_data[debut-1])
		liste_Z, liste_N, liste_E, liste_dE, liste_R0, liste_sig, liste_Rch, liste_Rp = [],[],[],[],[],[],[],[]
		for ligne in self.liste_data[debut:debut+N]:
			brokenLine = re.split("\n",ligne)
			brokenLine = re.split(self.separator,brokenLine[0])
			strippedLine = filter(not_empty, brokenLine)
			liste_Z.append(strippedLine[0])
			liste_N.append(strippedLine[1])
			liste_E.append(strippedLine[2])
			liste_dE.append(strippedLine[3])
			liste_R0.append(strippedLine[4])
			liste_sig.append(strippedLine[5])
			liste_Rch.append(strippedLine[6])
			liste_Rp.append(strippedLine[7])

		return { "Z": liste_Z, "N": liste_N, "E": liste_E, "dE": liste_dE, "R0": liste_R0, \
		         "sigma": liste_sig, "Rch": liste_Rch, "Rp": liste_Rp}


	# Obtain data for deformed nuclei
	def get_deformed(self):
		chaine = "DEFORMED ! Z   N     B [MeV]      dB [MeV]    b2 (SLy4) Meas."
		position = [ i for i, x in enumerate(self.liste_data) if x.find(chaine) > -1 ]
		debut = int(position[0])+2
		N = int(self.liste_data[debut-1])
		liste_Z, liste_N, liste_E, liste_dE, liste_beta, liste_eval = [],[],[],[],[],[]
		for ligne in self.liste_data[debut:debut+N]:
			brokenLine = re.split(self.separator,ligne)
			strippedLine = filter(not_empty, brokenLine)
			liste_Z.append(strippedLine[0])
			liste_N.append(strippedLine[1])
			liste_E.append(strippedLine[2])
			liste_dE.append(strippedLine[3])
			liste_beta.append(strippedLine[4])
			liste_eval.append(strippedLine[5])

		return { "Z": liste_Z, "N": liste_N, "E": liste_E, "dE": liste_dE, "beta": liste_beta, \
		         "evaluation": liste_eval}

	# Obtain data for neutron odd-even mass differences
	def get_delta3n(self):
		chaine = "DELTA3_N !   Z   N      Delta3       Error  Meas."
		position = [ i for i, x in enumerate(self.liste_data) if x.find(chaine) > -1 ]
		debut = int(position[0])+2
		N = int(self.liste_data[debut-1])
		liste_Z, liste_N, liste_D3n, liste_err, liste_eval = [],[],[],[],[]
		for ligne in self.liste_data[debut:debut+N]:
			brokenLine = re.split(self.separator,ligne)
			strippedLine = filter(not_empty, brokenLine)
			liste_Z.append(strippedLine[0])
			liste_N.append(strippedLine[1])
			liste_D3n.append(strippedLine[2])
			liste_err.append(strippedLine[3])
			liste_eval.append(strippedLine[4])

		return { "Z": liste_Z, "N": liste_N, "D3n": liste_D3n, "err": liste_err, "evaluation": liste_eval}

	# Obtain data for proton odd-even mass differences
	def get_delta3p(self):
		chaine = "DELTA3_P !   N   Z      Delta3       Error  Meas."
		position = [ i for i, x in enumerate(self.liste_data) if x.find(chaine) > -1 ]
		debut = int(position[0])+2
		N = int(self.liste_data[debut-1])
		liste_Z, liste_N, liste_D3p, liste_err, liste_eval = [],[],[],[],[]
		for ligne in self.liste_data[debut:debut+N]:
			brokenLine = re.split(self.separator,ligne)
			strippedLine = filter(not_empty, brokenLine)
			liste_N.append(strippedLine[0])
			liste_Z.append(strippedLine[1])
			liste_D3p.append(strippedLine[2])
			liste_err.append(strippedLine[3])
			liste_eval.append(strippedLine[4])

		return { "Z": liste_Z, "N": liste_N, "D3p": liste_D3p, "err": liste_err, "evaluation": liste_eval}

	# Obtain data for proton odd-even mass differences
	def get_isomers(self):
		chaine = "SDSTATES !   Z   N     B [MeV]   ESD [MeV]  beta   Reference"
		position = [ i for i, x in enumerate(self.liste_data) if x.find(chaine) > -1 ]
		debut = int(position[0])+2
		N = int(self.liste_data[debut-1])
		liste_Z, liste_N, liste_E, liste_ESD, liste_b2 = [],[],[],[],[]
		for ligne in self.liste_data[debut:debut+N]:
			brokenLine = re.split(self.separator,ligne)
			strippedLine = filter(not_empty, brokenLine)
			liste_Z.append(strippedLine[0])
			liste_N.append(strippedLine[1])
			liste_E.append(strippedLine[2])
			liste_ESD.append(strippedLine[3])
			liste_b2.append(strippedLine[4])

		return { "Z": liste_Z, "N": liste_N, "E": liste_E, "ESD": liste_ESD, "beta": liste_b2}
	
	# Obtain data for deformed nuclei
	def get_meas2013(self):
		chaine = "MEAS2013 ! Z   N     B [MeV]      dB [MeV]    Q2 (UNE1) Meas."
		position = [ i for i, x in enumerate(self.liste_data) if x.find(chaine) > -1 ]
		debut = int(position[0])+2
		N = int(self.liste_data[debut-1])
		liste_Z, liste_N, liste_E, liste_dE, liste_beta, liste_eval = [],[],[],[],[],[]
		for ligne in self.liste_data[debut:debut+N]:
			brokenLine = re.split(self.separator,ligne)
			strippedLine = filter(not_empty, brokenLine)
			#print "DEBUG: ", strippedLine
			liste_Z.append(strippedLine[0])
			liste_N.append(strippedLine[1])
			liste_E.append(strippedLine[2])
			liste_dE.append(strippedLine[3])
			liste_beta.append(strippedLine[4])
			liste_eval.append(strippedLine[5])

		return { "Z": liste_Z, "N": liste_N, "E": liste_E, "dE": liste_dE, "beta": liste_beta, \
		         "evaluation": liste_eval}

