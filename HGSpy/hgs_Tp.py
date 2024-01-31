'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS Adiabatic flame temperature algorithm

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import numpy as np

from .hgs        import HGSData
from .cr         import cr
from .utils      import raiseError
from .hgs_prop   import hgs_prop_ids
from .hgs_eq     import options as opt_eq
from .hgs_solver import hgs_solver, options as opt_sec


def hgs_Tp_ids(ids, n0, typ, V0, P, flow, solver, Tstar, opt_eq, opt_sci, opt_sec, hgs_data):
	'''
	Main function for hgs_Tp working with ids instead of species
	'''
	# Compute initial enthalpy
	H = hgs_prop_ids(ids,n0,V0,P,['H'],hgs_data)[0] if typ == 'T' else V0[0]
	# Solve for adiabatic flame temperature
	Tp, n, flag = hgs_solver(ids,n0,'H',H,P,solver=solver,Tstar=Tstar,
		opt_sec=opt_sec,opt_eq=opt_eq,hgs_data=hgs_data)
	return Tp, n, flag

@cr('HGS.Tp')
def hgs_Tp(species, n0, typ, V0, P, flow='shifting', solver='hgs_secant', Tstar=3000, 
	opt_eq=opt_eq, opt_sci={}, opt_sec=opt_sec, hgs_data=HGSData.load()):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	Tp, n, species, flag = hgs_tp(species, n0, tipo, V0, P, **kwargs)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_tp calculates the reaction temperature considering dissociation,
	and the products composition in equilibrium

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	species --> String or numbers of species
	n0 --> [mols] Number of mols of each species
	type --> Entry type that defines the state of the input.
			 It can be 'T' or 'H'
	V0 --> Entry that should be for type:'T'   V0=T [K] input temperature
										 'H'   V0=H [kJ] input enthalpy
	P --> [bar] Mixture pressure
	**kwargs --> opt_eq = Options for the minimize Scipy function
				 opt_sec= Dictionary with the options for the secant method.
						"xmin": [K] Temperature minimum for the solver;
						"xmax" [K] Temperature maximum for the solver;
						"maxiter" Max iterations for the solver;
						"epsx" Diferential T where the solver reachs the solution;
						"epsy" Diferential S where the solver reachs the solution;
						"fchange" T difference where secant method is
								 changed by bisection method;
						"tipo" Select between: 'Frozen' for frozen flow
											  'Shifting' for shifting flow
						"info" Detailed info == 1; No info == 0.
						"dTp" Improve the velocity with the approximation of
							  parabola. +- dTp
						opt_sec = {"xmin": 300, "xmax": 4000, "maxiter": 200,
								   "epsx": 0.1, "epsy": 1, "tipo": "Shifting",
								   "fchange": 5, "info": 0, "dTp": 100}

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	Tp --> [K] Final temperature
	n --> [mol] Final mixture
	species --> String or numbers of species
	flag --> Solver error detection:
				  1  Solver has reached the solution
				 -1  Solver failed. Maximum iterations
				 -2  Solver failed. Initial sign change not found

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	if typ not in ['H','T']:                         raiseError(f'Wrong type = {type}')
	if type(V0) in (float,int,np.float64,np.float32): V0 = [V0]*len(species)
	if len(V0) == 1:                                  V0 = [V0[0]]*len(species)

	ids = hgs_data.id(species)
	# Rebuild mixtures
	if np.max(ids) >= len(hgs_data):
		species, n0, V0 = hgs_data.rebuild(species,n0,V0)
		ids             = hgs_data.id(species)
		
	Tp, n, flag = hgs_Tp_ids(ids, n0, typ, V0, P, flow, solver, Tstar, opt_eq, opt_sci, opt_sec, hgs_data)

	# Return output
	return Tp, n, species, flag