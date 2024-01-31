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
from .utils      import raiseError, raiseWarning
from .hgs_prop   import hgs_prop_ids
from .hgs_eq     import hgs_eq_ids, options as opt_eq
from .hgs_solver import hgs_solver, options as opt_sec
from .hgs_secant import hgs_secant


def hastobeS_shifting(S,Tstar,Pstar,ids,ni,opt_eq,hgs_data):
	'''
	Function to look for entropy given a pressure and temperature
	'''
	ni,_ = hgs_eq_ids(ids,ni,[Tstar]*len(ni),Pstar,opt_eq,hgs_data)
	Si = hgs_prop_ids(ids,ni,[Tstar]*len(ni),Pstar,['S'],hgs_data)[0]
	return Si - S, ni

def hastobeM_shifting(V1,m,h1,Pstar,S,ids,n,opt_eq,hgs_data):
	'''
	Function to look for Mach number given a pressure
	'''
	Tstar, n, flag = hgs_secant(lambda Ti, ni: hastobeS_shifting(S,Ti,Pstar,ids,ni,opt_eq,hgs_data),n,opt_sec,ilevel=1)
	if flag != 1: raiseError('uhhh error in hastobeM HGSsecant failed flag=%d'%flag)
	a, H2 = hgs_prop_ids(ids,n,[Tstar]*len(n),Pstar,['a','H'],hgs_data)
	v2    = np.sqrt(2*1000*(h1-H2/m))
	M1    = v2/a
	#print('Pstar=%e Tstar=%e n(3)=%e V1(goal)=%e M=%e \n'%(Pstar,Tstar,n(3),V1,M1))
	return V1 - M1, n

def hastobeS_frozen(S,Tstar,Pstar,ids,ni,opt_eq,hgs_data):
	'''
	Function to look for entropy given a pressure and temperature
	'''
	Si = hgs_prop_ids(ids,ni,[Tstar]*len(ni),Pstar,['S'],hgs_data)[0]
	return Si - S, ni

def hastobeM_frozen(V1,m,h1,Pstar,S,ids,n,opt_eq,hgs_data):
	'''
	Function to look for Mach number given a pressure
	'''
	Tstar, n, flag = hgs_secant(lambda Ti, ni: hastobeS_frozen(S,Ti,Pstar,ids,ni,opt_eq,hgs_data),n,opt_sec,ilevel=1)
	if flag != 1: raiseError('uhhh error in hastobeM HGSsecant failed flag=%d'%flag)
	a, H2 = hgs_prop_ids(ids,n,[Tstar]*len(n),Pstar,['a','H'],hgs_data)
	v2    = np.sqrt(2*1000*(h1-H2/m))
	M1    = v2/a
	#print('Pstar=%e Tstar=%e n(3)=%e V1(goal)=%e M=%e \n'%(Pstar,Tstar,n(3),V1,M1))
	return V1 - M1, n


def hgs_isentropic_ids(ids, n0, T0, P0, type, V1, flow, solver, Tstar, 
	opt_eq, opt_sci, opt_sec, hgs_data):
	'''
	Main function for hgs_isentropic working with ids instead of species
	'''
	# Compute initial entropy and enthalpy
	S, Mm1, H1 = hgs_prop_ids(ids,n0,T0,P0,['S','Mm','H'],hgs_data) # Inlet  properties
	m1         = np.sum(n0)*Mm1*1e-3
	h1         = H1/m1

	if type == 'P':
		P1 = V1
		Tp, n, flag = hgs_solver(ids,n0,'S',S,P1,
			flow     = flow,
			solver   = solver,
			Tstar    = Tstar,
			opt_sec  = opt_sec,
			opt_eq   = opt_eq,
			hgs_data = hgs_data
		)
		Mm2, a2, H2 = hgs_prop_ids(ids,n,[Tp]*len(ids),P1,['Mm','a','H'],hgs_data) # Outlet properties
		m2          = np.sum(n)*Mm2*1e-3
		h2          = H2/m2
		v2          = np.sqrt(2*1000*(h1 - h2)) # Enthalpy must be en J / kg !
		V2          = v2/a2
	elif type == 'M':
		M2 = V1
		# To avoid problems with imaginary numbers
		hastobeM = hastobeM_shifting
		if flow.lower() == 'frozen': 
			P0 *= 0.9
			hastobeM = hastobeM_frozen
		# Update settings to find pressure change
		opt_sec2 = opt_sec.copy()
		opt_sec2['xmin']    = 0.1
		opt_sec2['xmax']    = P0
		opt_sec2['maxiter'] = 50
		opt_sec2['epsx']    = 0.01
		opt_sec2['epsy']    = 0.001
		opt_sec2['fchange'] = 1
		# Pressure search, imposing M
		P1,n,flag = hgs_secant(lambda Pi, ni : hastobeM(V1,m1,h1,Pi,S,ids,ni,opt_eq,hgs_data),n0,opt_sec2)
		# Flag error return
		if flag != 1:
			raiseWarning('HGSisentropic error in finding pressure')
		# T calculation
		Tp, n, flag = hgs_solver(ids,n,'S',S,P1,
			flow     = flow,
			solver   = solver,
			Tstar    = Tstar,
			opt_sec  = opt_sec,
			opt_eq   = opt_eq,
			hgs_data = hgs_data
		)
		a2 = hgs_prop_ids(ids,n,[Tp]*len(ids),P1,['a'],hgs_data)[0]  # Outlet properties
		v2 = M2*a2
		V2 = P1

	return Tp, n, v2, V2, flag

@cr('HGS.isentropic')
def hgs_isentropic(species, n0, T0, P0, typ, V1, flow='shifting', solver='hgs_secant', Tstar=3000,
	opt_eq=opt_eq, opt_sci={}, opt_sec=opt_sec, hgs_data=HGSData.load()):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	Tp, n, species, v2, M2, flag = hgs_isentropic(species, n0, T0, P0, P1, **kwargs)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_tp calculates the outlet variables for an isentropic expansion

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-------------------------------------------------------------------------------
	species --> String or code of species
	n0 --> [mols] Number of mols of each species
	T0 --> [K] Initial temperature
	P0 --> [bar] Inlet pressure
	type --> Entry type that defines the state of the input. 
				  It can be 'P' or 'M'
	V1 --> Value for type:'P'   V1=P [bar] output pressure
						  'M'   V1=M [] output Mach. Has to be >=1
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
	Tp --> [K] Exit temperature
	n --> [mols] Species resultant mols
	species --> String or numbers of species
	v2 --> [m/s] Velocity of the mixture
	M2 --> [Mach] Mach of the mixture
	flag --> Solver error detection:
				  1  Solver has reached the solution
				 -1  Solver failed. Maximum iterations
				 -2  Solver failed. Initial sign change not found

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	if typ not in ['P','M']:                          raiseError(f'Wrong type = {typ}')
	if type(T0) in (float,int,np.float64,np.float32): T0 = [T0]*len(species)
	if len(T0) == 1:                                  T0 = [T0[0]]*len(species)

	ids = hgs_data.id(species)
	# Rebuild mixtures
	if np.max(ids) >= len(hgs_data):
		species, n0, T0 = hgs_data.rebuild(species,n0,T0)
	
	Tp, n, v2, V2, flag = hgs_isentropic_ids(ids,n0,T0,P0,typ,V1,flow,solver,Tstar,opt_eq,opt_sci,opt_sec,hgs_data)

	return Tp, n, species, v2, V2, flag