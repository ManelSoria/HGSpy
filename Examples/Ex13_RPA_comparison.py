#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************%
#
# RPA vs HGStp && HGSisentropic.
# LH2-LOX reaction
#
# RPA is a code for rocket combustion and expansion in a nozzle.
# Search it on the web for more info. https://www.rocket-propulsion.com/index.htm
# We compare our HGS algorithms with their results for a LH2-LOX combustion
#
# You can find a screen capture of the RPA solution in the file
# RPA_H2_O2_ex.PNG
from __future__ import print_function

import numpy as np, HGSpy as HGS


HGS.set_options('warnings',False) # Deactivate warnings

# RPA data (ignore H2O2 and HO2 species)
# Combustion data
RPATcomb = 3004.0779
RPAcomb  = [0.4730081, 0.0000519, 0.5040833, 0.0063140, 0.0001233, 0.0164185]
# Exit of the isentropic expansion
RPATis   = 1534.1342
RPAis    = [0.4833980, 0, 0.4833980, 0.0000007, 0, 0.0000183]


# Basic info
of_ratio = 4.1 # The OF Ratio is done with kg oxidizer/ kg fuel
species  = ['H2', 'O2', 'H2O', 'OH', 'O', 'H']
n        = np.zeros((len(species),))
Tin      = 90  # [K]
Tref     = 300 # [K]
Pin      = 45  # [bar]
Pexit    = 1   # [bar]

mH2      = 1   # [kg]
mO2      = mH2*of_ratio # [kg]

# Convert kg to mols
MmH2 = HGS.single('H2','Mm',1,1)
MmO2 = HGS.single('O2','Mm',1,1)

n[0] = mH2/(MmH2/1000)
n[1] = mO2/(MmO2/1000)

# First of all, we need to calculate de inlet enthalpy using INIST
# (only the values) from a reference T of 300 K.
# deltaH_H2 = (INIST('H2','h_pt',45,300) - INIST('H2','h_pt',45,90)) * mH2
# deltaH_O2 = (INIST('O2','h_pt',45,300) - INIST('O2','h_pt',45,90)) * mO2

deltaH_H2 = 2.8577525044e+03*mH2
deltaH_O2 = 393.3914211*mO2

H_H2_HGS  = HGS.single('H2','h',Tref,Pin)*n[0]
H_O2_HGS  = HGS.single('O2','h',Tref,Pin)*n[1]

# Inlet enthalpy
Hin = H_H2_HGS - deltaH_H2 + H_O2_HGS - deltaH_O2

# Now, we run the combustion with HGStp
Tcomb,ncomb,_,_  = HGS.Tp(species,n,'H',Hin,Pin)

# Now, we run the isentropic expansion
Tis,nis,_,M2,_,_ = HGS.isentropic(species,ncomb,Tcomb,Pin,'P',Pexit)

# Comparing results
print('After the combustion')
print('RPA\t %.3f K\t\t HGS: %.3f K'%(RPATcomb,Tcomb))
print('Molar fractions')
for ii,s in enumerate(species):
	print('%s\t RPA:%.3e \t\t HGS: %.3e'%(s,ncomb[ii]/np.sum(ncomb),RPAcomb[ii]))

print()
print('After the expansion')
print('RPA\t %.3f K\t\t HGS: %.3f K'%(RPATis,Tis))
print('Molar fractions')
for ii,s in enumerate(species):
	print('%s\t RPA:%.3e \t\t HGS: %.3e'%(s,nis[ii]/np.sum(nis),RPAis[ii]))

HGS.cr_info()