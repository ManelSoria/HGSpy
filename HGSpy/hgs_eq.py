'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS Equilibrium algorithm

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import numpy as np
from scipy.optimize import minimize, Bounds

from .cr       import cr
from .hgs      import HGSData
from .utils    import raiseError, raiseWarning
from .hgs_prop import hgs_prop_ids


options = {
	'method':'SLSQP',
	'tol':1e-9,
	'options':{'disp':False},
}


# Parameters minimization
def parameters_min(ids, n0, hgs_data):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	bounds, linear = parameters_min(ids, n0)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	parameters_min returns the parameters for minimize "SQLP" solver

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	id --> Id of species
	n0 --> [mol] Species mols

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	bounds --> Bound limits
	linear --> Linear equality constrains

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	# Bounds
	bounds = Bounds([0]*len(ids),[np.inf]*len(ids))

	# Equality
	Elems = np.unique([e for i in ids for e in hgs_data['ena'][i]])

	Aeq = np.zeros((len(Elems),len(ids)),np.double)
	beq = np.zeros((len(Elems),),np.double)

	for ii,e in enumerate(Elems):
		for jj,q in enumerate(ids):
			for kk,ee in enumerate(hgs_data['ena'][q]):
				if ee == e:
					beq[ii]    += n0[jj]*hgs_data['nat'][q][kk]
					Aeq[ii][jj] = hgs_data['nat'][q][kk]

	linear = {"type": "eq","fun": lambda x: np.dot(Aeq,x) - beq}
	return bounds, linear


def hgs_eq_ids(ids, n0, T, P, options, hgs_data):
	'''
	Main function for hgs_eq working with ids instead of species
	'''
	# Function to minimize
	minG = lambda x: hgs_prop_ids(ids,x,T,P,['g'],hgs_data)[0]

	# Function minimization parameter
	bounds, linear = parameters_min(ids,n0,hgs_data)
	res = minimize(minG,n0,method=options['method'],tol=options['tol'],
				   constraints=linear,options=options['options'],bounds=bounds)

	if not res.success:
		raiseWarning("Ups,... minimize has failed in hgs_eq.")

	return res.x, res.fun

@cr('HGS.eq')
def hgs_eq(species, n0, T, P, options=options, hgs_data=HGSData.load()):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	species, n, Gmin = hgs_eq(species, n0, T, P, **kwargs)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_eq calculates the species mols equilibrium at a certain temperature
	and pressure

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	species --> String or numbers of species
	n0 --> [mol] Initial mixture
	T --> [K] Temperature. Could be a single value or an array.
	P --> [bar] Pressure
	**kwargs --> opti_eq= Options for the minimize Scipy function

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	species --> Species
	n --> [mol] Final mixture
	Gmin --> [kJ] Minimum Gibbs free energy

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	if len(species) != len(n0):
		raiseError("Ups..., Species and mols have not the same length. Check it")
	if type(T) in (float,int,np.float64,np.float32): T = [T]*len(species)
	if len(T) == 1:                                  T = [T[0]]*len(species)

	ids = hgs_data.id(species)
	# Rebuild mixtures
	if np.max(ids) >= len(hgs_data):
		species, n0, _ = hgs_data.rebuild(species,n0,T)

	return species, *hgs_eq_ids(ids, n0, T, P, options, hgs_data)