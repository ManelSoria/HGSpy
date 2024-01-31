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

from .hgs            import HGSData
from .cr             import cr
from .utils          import raiseError, raiseWarning
from .definitions    import g0
from .hgs_prop       import hgs_prop_ids
from .hgs_eq         import options as opt_eq
from .hgs_solver     import hgs_solver, options as opt_sec
from .hgs_isentropic import hgs_isentropic_ids


@cr('HGS.nozzle')
def hgs_nozzle(species, n0, T0, P0, P, Pa, flow='shifting', solver='hgs_secant', Tstar=3000, 
	opt_eq=opt_eq, opt_sci={}, opt_sec=opt_sec, hgs_data=HGSData.load()):
	'''
	**************************************************************************
	
	 [species,n,T,v,M,A,F,Isp] = HGSnozzle(species,n0,T0,P0,Pe,Pa,Fro_Shift,
								 options1,options2)
	
	**************************************************************************
	 
	 HGSnozzle evaluates different flow properties as a function of the
	 pressure, during a isentropic expansion beginning with a very low
	 velocity
	  
	**************************************************************************
	 Inputs:
	--------------------------------------------------------------------------
	 species --> String or code of inlet species
	 n0 --> [mols] Number of mols/s of each inlet species
	 T0 --> [K] Inlet temperature
	 P0 --> [bar] Inlet pressure
	 P -->  [bar] Pressure vector (with all the pressures that have to be
			 evaluated)
	 Pa --> [bar] Atmospheric pressure
	 Fro_Shift --> Select between: 'Frozen' for frozen flow
								   'Shifting' for shifting flow
	
	 Outputs:
	--------------------------------------------------------------------------
	 species --> String or code of species
	 n --> [mols] Matrix of pecies mols, sorted as: n(species, pressure)
	 T --> [K] Exit temperature
	 M --> [] Exit Mch
	 A --> [m^2] Exit area
	 F --> [N] Thrust
	 Isp --> [s^]Specific impulse, g0 = 9.807 m/s^2
	
	**************************************************************************
	* Python HGS 1.0 from Matlab HGS 2.1
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	'''
	if type(T0) in (float,int,np.float64,np.float32): T0 = [T0]*len(species)
	if len(T0) == 1:                                  T0 = [T0[0]]*len(species)

	ids = hgs_data.id(species)
	# Rebuild mixtures
	if np.max(ids) >= len(hgs_data):
		species, n0, T0 = hgs_data.rebuild(species,n0,T0)

	# Total mass
	mm  = hgs_prop_ids(ids,n0,T0,P0,['Mm'],hgs_data)[0] # g/mol
	m   = np.sum(n0)*mm*1e-3 # kg/s

	# Preallocate
	P   = np.array(P)
	T   = np.zeros_like(P)
	n   = np.zeros((len(species),len(P)))
	M   = np.zeros_like(P)
	v   = np.zeros_like(P)
	A   = np.zeros_like(P)
	F   = np.zeros_like(P)
	Isp = np.zeros_like(P)

	# Run loop
	for ii in range(len(P)):
		print('P = %f,  %i /%i'%(P[ii],ii+1,len(P)))
		T[ii],n[:,ii],_,M[ii],flag = hgs_isentropic_ids(ids,n0,T0,P0,'P',P[ii],flow,solver,Tstar,opt_eq,opt_sci,opt_sec,hgs_data)
		if not flag == 1: raiseWarning('HGSnozzle failed to converge/1 flag=%d',flag) 
		Rg,a,Mm = hgs_prop_ids(ids,n[:,ii],[T[ii]]*len(ids),P[ii],['Rg','a','Mm'],hgs_data) # kJ/(kg*K), m/s, g/mol
		rho     = P[ii]*1e5/(Rg*1000*T[ii])        # kg/m^3 Convert bar to Pa g 2 kG
		v[ii]   = M[ii]*a                          # m/s
		A[ii]   = m/(v[ii]*rho)                    # m^2
		F[ii]   = m*v[ii] + A[ii]*(P[ii] - Pa)*1e5 # Convert bar to Pa
		Isp[ii] = v[ii]/g0

	return species, n, T, v, M, A, F, Isp