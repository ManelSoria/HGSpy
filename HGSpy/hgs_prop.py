'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS Prop main functions

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import numpy as np

from .hgs         import HGSData
from .cr          import cr
from .utils       import raiseError
from .definitions import R


## ---------- Properties ---------- #
def rg(mm):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	Rg = prop_rg(mm)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	prop_rg calculates the specific R from a mixture

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	mm --> [g/mol] Molar mass

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	Rg --> [kJ/(kg*K)] Specific R

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	return R/mm*1000.

def gamma(cp,cv):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	gamma = prop_gamma(cp, cv)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	prop_gamma calculates adiabatic expansion coefficient from cp and cv

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	cp --> [kJ/K] Constant pressure coefficient
	cv --> [kJ/K] Constant volume coefficient

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	gamma --> adiabatic expansion coefficient

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	return cp/cv

def sound(gamma,rg,T):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	a = prop_a(gamma, rg, T)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	prop_gamma calculates the sound velocity

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	gamma --> adiabatic expansion coefficient
	Rg --> [kJ/(kg*K)] Specific R
	T --> [K] Temperature

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	a --> [m/s] Sound speed

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	return np.sqrt(gamma*T*rg*1.e3)

def partial(n, P, ids, hgs_data):
	'''
	Compute partial pressure of gas mixtures
	'''
	# Check that  the mixture only contains gases
	st = [hgs_data['state'][i] == 'G' for i in ids]
	if not np.all(st):
		raiseError('Ups,.. Right now entropy can be calculated only for gas mixtures')
	# Compute partial pressure
	nt = np.sum(n)
	return [P*nn/nt for nn in n]


## ---------- Function   ---------- #
def hgs_prop_ids(ids, n, T, P, args, hgs_data):
	'''
	Main function for hgs_prop working with ids instead of species
	'''
	Tm  = np.dot(T,n)/np.sum(n) # Average temperature

	# args: (str) - Property(need to be calculated before)
	#       (mm)  - Molar mass()
	#       (cp)  - Cp(10 - Burcat)
	#       (cv)  - Cv(10 - Burcat & 2 - Cp)
	#        (h)  - H(10 - Burcat)
	#        (s)  - S(10 - Burcat)
	#        (g)  - G(10 - Burcat)
	#       (Rg)  - Rg(1 - Mm)
	#    (gamma)  - Gamma(2 - Cp & 3 - Cv)
	#        (a)  - Sound Velocity(1 - Mm & 2 - Cp & 3 - Cv & 7 - Rg & 8 - Gamma)
	#     (coef)  - Burcat Coef()
	if len(args) == 0:
		args = ['mm','cp','cv','h','s','g','rg','gamma','a']

	# Generate output
	out = []
	for a in args:
		if a.lower() == 'mm': # Molar mass
			mm = [hgs_data['mm'][i] for i in ids]
			nt = np.sum(n)
			out.append(np.dot(n,mm)/nt)
		if a.lower() == 'cp': # Cp
			cp = [hgs_data.cp(i,Ti) for i,Ti in zip(ids,T)]
			out.append(np.dot(n,cp))
		if a.lower() == 'cv': # Cv
			cv = [hgs_data.cv(i,Ti) for i,Ti in zip(ids,T)]
			out.append(np.dot(n,cv))
		if a.lower() == 'h': # H
			h = [hgs_data.h(i,Ti) for i,Ti in zip(ids,T)]
			out.append(np.dot(n,h))
		if a.lower() == 's': # S
			P_i = partial(n, P, ids, hgs_data)
			s   = [hgs_data.s(i,Ti,p) for i,p,Ti in zip(ids,P_i,T)]
			out.append(np.dot(n,s))
		if a.lower() == 'g': # G
			P_i = partial(n, P, ids, hgs_data)
			g   = [hgs_data.g(i,Ti,p) for i,p,Ti in zip(ids,P_i,T)]
			out.append(np.dot(n,g))
		if a.lower() == 'rg': # Rg
			mm = [hgs_data['mm'][i] for i in ids]
			nt = np.sum(n)
			out.append(rg(np.dot(n,mm)/nt))
		if a.lower() == 'gamma': # Gamma
			cp = [hgs_data.cp(i,Ti) for i,Ti in zip(ids,T)]
			cv = [hgs_data.cv(i,Ti) for i,Ti in zip(ids,T)]
			out.append(gamma(np.dot(n,cp),np.dot(n,cv)))
		if a.lower() == 'a': # a (Sound velocity)
			mm = [hgs_data['mm'][i] for i in ids]
			nt = np.sum(n)
			cp = [hgs_data.cp(i,Ti) for i,Ti in zip(ids,T)]
			cv = [hgs_data.cv(i,Ti) for i,Ti in zip(ids,T)]
			out.append(sound(gamma(np.dot(n,cp),np.dot(n,cv)),rg(np.dot(n,mm)/nt),Tm))
		if a.lower() == 'coef': # Burcat coefficients
			out.append([hgs_data.coef(i,Ti) for i,Ti in zip(ids,T)])

	return out

@cr('HGS.single')
def hgs_single(species, prop, T, P, hgs_data=HGSData.load()):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	res = hgs_single(species, prop, T, P)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_single returns the property of a species

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	species --> String or numbers of species
	prop --> Property requested (see below)
	T --> [K] Temperature
	P --> [bar] Pressure

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	res --> Property result
		  mm [g/mol]
		  cp [kJ/(mol*K)]
		  cv [kJ/(mol*K)]
		  h [kJ/mol]
		  s [kJ/(mol*K)]
		  g [kJ/mol]

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	if not prop.lower() in ['mm','cp','cv','h','s','g','rg','gamma','a','coef']:
		raiseError(f'Property {prop} not understood!')
	if type(T) in (float,int,np.float64,np.float32): T = [T]

	ids = hgs_data.id(species)
	n   = [1]
	# Rebuild mixtures
	if np.max(ids) >= len(hgs_data):
		species, n, T = hgs_data.rebuild(species,[1],T)
		ids           = hgs_data.id(species)
	
	return hgs_prop_ids(ids, n, T, P, [prop], hgs_data)[0]

@cr('HGS.prop')
def hgs_prop(species, n, T, P, *args, hgs_data=HGSData.load()):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	var = hgs_prop(species, n, T, P, *args)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_prop returns the properties of the mixture of gasses

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	species --> String or numbers of species
	prop --> Property requested (see below)
	T --> [K] Temperature
	P --> [bar] Pressure
	*args --> Expected return: 'Mm' 'Cp' 'Cv' 'H' 'S' 'G' 'Rg' 'gamma' 'a'
							   If it is empty, all the properties will be
							   return

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	var --> Property result
		  mm [g/mol]
		  cp [kJ/K]
		  cv [kJ/K]
		  h [kJ]
		  s [kJ/K]
		  g [kJ]
		  Rg [kJ/(kg*K)]
		  gamma
		  a [m/s]

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	if type(species) is str:                         species = [species]
	if type(n) in (float,int,np.float64,np.float32): n       = [n]
	if type(T) in (float,int,np.float64,np.float32): T       = [T]*len(species)
	if len(T) == 1:                                  T       = [T[0]]*len(species)

	if len(n) is not len(species):
		raiseError('Ups..., Species (%d) and mols (%d) lengths are not the same!'%(len(species),len(n)))

	ids = hgs_data.id(species)
	# Rebuild mixtures
	if np.max(ids) >= len(hgs_data):
		species, n, T = hgs_data.rebuild(species,n,T)
		ids           = hgs_data.id(species)

	# Return properties
	return hgs_prop_ids(ids, n, T, P, args, hgs_data)