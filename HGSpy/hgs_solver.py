'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS Solvers

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import numpy as np, scipy.optimize

from .hgs        import HGSData
from .cr         import cr, cr_start, cr_stop
from .utils      import raiseError
from .hgs_secant import hgs_secant
from .hgs_eq     import hgs_eq_ids, options as opt_eq
from .hgs_prop   import hgs_prop_ids


options = {
	"xmin":    300,
	"xmax":    4000,
	"maxiter": 200, 
	"epsx":    5, 
	"epsy":    1, 
	"fchange": 500, 
	"info":    0,
	"dTp":     100,
}


def hastobezeroH_shifting(T, P, ni, ids, V0, opt_eq, hgs_data):
	'''
	Function that must be zero for tipo == H
	'''
	ni,_ = hgs_eq_ids(ids,ni,[T]*len(ids),P,opt_eq,hgs_data)
	return hgs_prop_ids(ids,ni,[T]*len(ids),P,['H'],hgs_data)[0] - V0, ni

def hastobezeroH_frozen(T, P, ni, ids, V0, opt_eq, hgs_data):
	'''
	Function that must be zero for tipo == H
	'''
	return hgs_prop_ids(ids,ni,[T]*len(ids),P,['H'],hgs_data)[0] - V0, ni

def hastobezeroS_shifting(T, P, ni, ids, V0, opt_eq, hgs_data):
	'''
	Function that must be zero for tipo == S
	'''
	ni,_ = hgs_eq_ids(ids,ni,[T]*len(ids),P,opt_eq,hgs_data)
	return hgs_prop_ids(ids,ni,[T]*len(ids),P,['S'],hgs_data)[0] - V0, ni

def hastobezeroS_frozen(T, P, ni, ids, V0, opt_eq, hgs_data):
	'''
	Function that must be zero for tipo == S
	'''
	return hgs_prop_ids(ids,ni,[T]*len(ids),P,['S'],hgs_data)[0] - V0, ni


@cr('HGS.solver')
def hgs_solver(ids, n0, typ, V0, P, flow='shifting', solver='hgs_secant', Tstar=3000, 
	opt_eq=opt_eq, opt_sci={}, opt_sec=options, hgs_data=HGSData.load()):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	Tp, n, species, flag = hgs_eqcond(species, n0, tipo, V0, P, **kwargs)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_secant computes temperature and reaction products that satisfy one of
	the following conditions
	if tipo=='H', the temperature Tp satisfies the condition H(Tp) = V0
	if tipo=='S', the temperature Tp satisfies the condition S(Tp) = V0

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	species --> Species name, ids are also accepted as entry
	n0 --> [mols] Number of mols of each species
	tipo --> Entry type. It coould be 'H' or 'S'
	V0 --> Entry that should be for tipo:'H'   V0=H [kJ]
										 'S'   V0=S [kJ/K]
	P --> [bar] Mixture pressure
	**kwargs --> opti_eq= Options for the minimize Scipy function
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
	if typ not in ['H','S']: raiseError(f'Wrong type = {typ}')

	if flow.lower() ==  'shifting':
		hastobezero = hastobezeroH_shifting if typ == 'H' else hastobezeroS_shifting
	else:
		hastobezero = hastobezeroH_frozen if typ == 'H' else hastobezeroS_frozen

	if solver == 'hgs_secant':
		Tp, n, flag = hgs_secant(lambda Ti, ni: hastobezero(Ti,P,ni,ids,V0,opt_eq,hgs_data),n0,opt_sec)
	else:
		# Use a solver from the scipy.optimize package
		Tp   = getattr(scipy.optimize,solver)(lambda Ti: hastobezero(Ti[0],P,n0,ids,V0,opt_eq,hgs_data)[0],Tstar,**opt_sci)
		n    = hastobezero(Tp[0],P,n0,ids,V0,opt_eq,hgs_data)[1]
		flag = 1
	return Tp, n, flag