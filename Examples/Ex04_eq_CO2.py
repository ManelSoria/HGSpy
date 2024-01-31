#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************
#
#  CO2 dissociation equilibrium solved with K
#
# 1 CO2    <-> 1 CO + (1/2)O2 
#   1-z          z       z/2
from __future__ import print_function

import numpy as np, HGSpy as HGS

from scipy.optimize import fsolve


Ru = HGS.R # kJ/molK
T  = 2500  # K
p  = 1     # bar

pref = 1   # bar (arbitrary pressure)

# Gibbs free energy (g) is computed for each species separately. The rest
# of properties are ignored.
gCO2 = HGS.prop('CO2',1,T,pref,'G')[0] # kJ
gCO  = HGS.prop('CO', 1,T,pref,'G')[0]
gO2  = HGS.prop('O2', 1,T,pref,'G')[0]

# Definition of K
deltag = gCO + 0.5*gO2 - gCO2
K      = np.exp(-deltag/(Ru*T))

# Amount of mols as a function of z
nCO  = lambda z: z
nO2  = lambda z: z/2
nCO2 = lambda z: 1-z

# Total amount as a function of z
nT = lambda z: nCO(z) + nO2(z) + nCO2(z)

# Molar fractions as a function of z
# Note: this code is probably very slow !
# It would be faster to use a single function for lhs(z,p)
xCO  = lambda z: nCO(z)/nT(z)
xO2  = lambda z: nO2(z)/nT(z)
xCO2 = lambda z: nCO2(z)/nT(z)

# From the definition of K the left hand side (lhs) of the equation is
# written
lhs = lambda z: (xCO(z)*xO2(z)**0.5 / xCO2(z) ) * (p/pref)**(1+1/2-1)

# The equation to solve is
eq = lambda z: lhs(z) - K
 
zi = 0.5 # a initial value (arbitrary) of z is given
z  = fsolve(eq,zi,xtol=1e-4,maxfev=4000)
print('K=%e z=%e'%(K,z))

print('x_CO =%6.4f'%xCO(z))
print('x_CO2=%6.4f'%xCO2(z))
print('x_O2 =%6.4f'%xO2(z))

# Print performance info
HGS.cr_info()