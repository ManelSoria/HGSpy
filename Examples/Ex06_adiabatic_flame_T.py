#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************
#
# Adiabatic flame temperature with dissociation
# this code is an example to understand HGStp
from __future__ import print_function

import numpy as np, HGSpy as HGS
import matplotlib.pyplot as plt

from scipy.optimize import fsolve


species = ['H2','O2','H2O','H','O','OH']
nr      = [   2,   1,    0,  0,  0,   0] # mol
Tr      = 350           # K 
P       = 10            # bar

Hin = HGS.prop(species,nr,Tr,P,'H')[0]


# Plot products enthalpy vs. T
T    = np.linspace(300,4000,20)
Hout = np.zeros_like(T)
for i in range(len(T)):
    _,comp,_ = HGS.eq(species,nr,T[i],P)
    Hout[i]  = HGS.prop(species,comp,T[i],P,'H')[0]

plt.figure()
plt.plot(T,Hout,'-or',label='Hout',linewidth=2)
plt.plot(T,Hin*np.ones_like(T),'-b',label='Hin',linewidth=2)
plt.legend() 
plt.xlabel('Temperature (K)')
plt.ylabel('Enthalpy (kJ/molK)')

# Function to be solved to find T so that deltaH=0
def deltaH(Tp):
    _,neq,_ = HGS.eq(species,nr,Tp,P) # find equilibrium composition
    Ho      = HGS.prop(species,neq,Tp,P,'H')[0]
    return Ho - Hin

# examples
print('deltaH @ 400 K = %.2f kJ/molK'%deltaH(400))
print('deltaH @ 4000 K = %.2f kJ/molK'%deltaH(4000))

# solving the equation
Tflame = fsolve(deltaH,3000)
print('Adiabatic flame temperature Tp = %.2f K'%Tflame)
HGS.cr_info() # Print performance info

plt.show()