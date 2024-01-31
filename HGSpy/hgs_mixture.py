'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS Mixture main functions

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import numpy as np

from .cr          import cr
from .utils       import raiseError


@cr('HGS.add_mixture')
def hgs_add_mixture(name, species, percent, hgs_data):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_add_mixture(name, species, percent)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_add_mixture add a mixture to the HGSdata database

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	name --> Name of the mixture f.e. 'Air'
	species --> species from the mixture f.e. {,...}
	percent --> Percentage of each species in the mixture [,...]

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	if not None in hgs_data.id(name,raise_error=False):
		raiseError('Ups,.. this name is already used')

	ids = hgs_data.id(species,raise_error=False)
	if None in ids:
		raiseError('Ups,.. at least one of the species is not in the hgs_data')

	if not np.sum(percent) == 100:
		raiseError('Ups,... The total percentage is not 100%')

	if type(species) is list and all(type(species[x]) is int for x in range(len(species))):
		spec = [""]*len(species)
		for ii in range(len(species)):
			spec[ii] = hgs_data["name"][ii]
		species = spec

	if name is not str:
		name = name[:]

	lm = len(hgs_data["mm"])
	if np.max(ids) > lm:
		buildS, buildP = [], []
		for ii in range(len(ids)):
			if ids[ii] > lm:
				for jj in range(len(hgs_data["cspec"][ids[ii]-lm])):
					buildS.append( hgs_data["cspec"][ids[ii]-lm][jj] )
					buildP.append( hgs_data["cper"][ids[ii]-lm][jj]*percent[ii]/100 )
			else:
				buildS.append( species[ii] )
				buildP.append( percent[ii] )

		species = buildS
		percent = buildP

	hgs_data["comb"].append(name)
	hgs_data["cspec"].append(species)
	hgs_data["cper"].append(percent)

	return hgs_data

@cr('HGS.subt_mixture')
def hgs_subt_mixture(name, hgs_data):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_subt_mixture(name)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_subt_mixture eliminates a mixture from the HGSdata database

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	name --> Name of the mixture f.e. 'Air'

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	ids = hgs_data.id(name,raise_error=False)
	if None in ids:
		raiseError('Ups,.. this name is not used in HGSdata')

	delet = ids[0] - len(hgs_data["mm"])
	hgs_data["comb"].pop(delet)
	hgs_data["cspec"].pop(delet)
	hgs_data["cper"].pop(delet)

	return hgs_data

@cr('HGS.rebuild')
def hgs_rebuild(species, n, T, hgs_data):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	species, n = hgs_rebuild(species, n)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_rebuild transform the mixtures in to their species. It also transforms
	the mols of the mixtures with the proportions assigned.

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	species --> String or numbers of species
	n --> [mols] Mols

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	species --> String of species
	n --> [mols] Mols

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	if type(species) is str: species = [species]
	if type(n) in (float,int,np.float64,np.float32): n = [n]
	if type(T) in (float,int,np.float64,np.float32): T = [T]*len(species)
	if len(T) == 1: T = [T[0]]*len(species)

	if type(species[0]) is int:
		nspec = [None]*len(species)
		for ii in range(len(species)):
			if species[ii] < len(hgs_data):
				nspec[ii] = hgs_data["name"][species[ii]]
			else:
				nspec[ii] = hgs_data["comb"][species[ii] - len(hgs_data)]
		species = nspec

	# Find the species of the mixture
	newspecies = []
	newn       = []
	newT       = []
	idstopop   = []
	for i,s in enumerate(species):
		ids = hgs_data.id(s)[0]
		if ids >= len(hgs_data):
			ids -= len(hgs_data)
			idstopop.append(i)
			newspecies += [spec for spec in hgs_data["cspec"][ids]]
			newn       += [n[i]*nsp/100 for nsp  in hgs_data["cper"][ids]]
			newT       += [T[i] for spec in hgs_data["cspec"][ids]]

	# Substitute mixtures by their species
	species  = [s  for i,s  in enumerate(species) if i not in idstopop]
	n        = [nn for i,nn in enumerate(n)       if i not in idstopop]
	T        = [tt for i,tt in enumerate(T)       if i not in idstopop]
	species += newspecies
	n       += newn
	T       += newT

	# Get the unique species
	Elem = np.unique(species)
	if not len(Elem) == len(species):
		nun = [0]*len(Elem)
		for ii,s in enumerate(species):
			for jj,e in enumerate(Elem):
				if e == s: nun[jj] += n[ii]
		n       = nun
		species = Elem.tolist()

	return species, n, T