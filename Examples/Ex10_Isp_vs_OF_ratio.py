#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************%
#
# Compute the ISP of H2, O2 mixture vs OF ratio (LIQUID FUEL)
from __future__ import print_function

import numpy as np, HGSpy as HGS
import matplotlib.pyplot as plt

from scipy.optimize import fsolve


HGS.set_options('warnings',False) # Deactivate warnings

species = ['H','H2','H2O','H2O2','HO2','O','O2','OH']

# Inlet temperature as if the reactives were gas at 300K
Te = 300 # K   reactives inlet temperature
Pc = 50  # bar chamber pressure 
P2 = 0.1 # bar nozzle exit    
	
vrof = [] # vector of OF ratios
visp = [] # vector of specific impulses

Tcs = 1800
T2s = 400

def DeltaH(Tc,ni_i):
	_,nc,_ = HGS.eq(species,ni_i,Tc,Pc)
	MMC,HC = HGS.prop(species,nc,Tc,Pc,'Mm','H')
	nC     = np.sum(nc)  # mixture total number of mols (1)
	mc     = nC*MMC*1e-3 # mixture mass kg
	hc     = HC/mc       # kJ/kgK
	return hc-h1

for rof in [2,3,4,5,6]:
	
	print('solving for rof=%f'%rof)

	# evaluate mol of each specie at inlet for a given ROF ratio
	nO2 = 1
	mO2 = nO2*32
	mH2 = mO2/rof
	nH2 = mH2/2
	
	ni_i = [
		0,   # H
		nH2,  # H2
		0,   # H2O
		0,   # H2O2
		0,   # HO2
		0,   # O
		nO2, # O2
		0    # OH
	];     

	ni_i = ni_i/np.sum(ni_i) # mole fractions  
	

	# Evaluate inlet properties with HGS assuming gas state
	MM,H = HGS.prop(species,ni_i,Te,Pc,'Mm','H')
	n    = np.sum(ni_i) # mixture total number of mols (1)
	m    = n*MM*1e-3    # mixture mass kg
	h1G  = H/m          # inlet mixture enthalpy in GAS state kJ/kgK 
						# we evaluate it just for comparision

	# Inlet enthalpy as if reactives were satured liquid at 10 bar
	# O2 (NIST) hv(404.36 K)-hl(119.62 K)=14.3753 kJ/mol
	# H2 (NIST) hv(413.96K)-hl(31.39K)=10.9495 kJ/mol

	# Enthalpy of O2 liq at Tsat 10 bar (kJ/mol)
	hO2 = HGS.single('O2','h',404.36,10) - 14.3753
	
	# Enthalpy of H2 liq at Tsat 10 bar (kJ/mol)
	hH2 = HGS.single('H2','h',413.96,10) - 10.9495 
	Hin = ni_i[1]*hH2 + ni_i[6]*hO2
	h1  = Hin/m # inlet mixture enthalpy in LIQUID state kJ/kg 

	# We find temperature at nozzle inlet solving for Delta_H=0
	# hgsTp function can't be used as it assumes gas state
	Tc  = fsolve(lambda T : DeltaH(T,ni_i),Tcs)
	print('Chamber outlet temperature Tc=%f'%Tc)
	Tcs = Tc # In next interation, we will begin with the solution obtained now
	
	_,ni_calc,_ = HGS.eq(species,ni_i,Tc,Pc)
	MM,S        = HGS.prop(species,ni_calc,Tc,Pc,'Mm','S')
	m           = np.sum(ni_calc)*MM*1e-3 # mixture mass kg (has to be as before!)
	s           = S/m

	# We use the previous value of T2 to begin the iterations
	# fzero is more robust solver in this case, but decreasing the
	# tolerance we can solve the problem with fsolve, that is faster
	T2,_,_,vt,_,_ = HGS.isentropic(species,ni_calc,Tc,Pc,'P',P2)
	T2s           = T2
	
	print('Nozzle outlet temperature T2=%f'%T2)
	
	Is = vt/9.81 # Is (optimal expansion, Pe=Pambient)    

	vrof.append(rof)
	visp.append(Is)
	
plt.plot(vrof,visp,linewidth=2)
plt.xlabel('OF ratio')
plt.ylabel('ISP (s)')
plt.title('ISP vs. OF ratio for a H2-O2 rocket')

HGS.cr_info()
plt.show()