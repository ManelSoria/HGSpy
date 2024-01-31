#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************
#
# H2O equilibrium dissociation for different values of
#             temperature (T)
#
# H20 <-> H2 + O2 + H + O + OH
from __future__ import print_function

import numpy as np, HGSpy as HGS
import matplotlib.pyplot as plt


n = 50
p = 1                       # bar
T = np.linspace(300,5000,n) # K

# loop to compute composition using hgseq
xH2  = np.zeros((n,),np.double)
xO2  = np.zeros((n,),np.double)
xH2O = np.zeros((n,),np.double)
xH   = np.zeros((n,),np.double)
xO   = np.zeros((n,),np.double)
xOH  = np.zeros((n,),np.double)
for i in range(n):
    print('Solving equilibrium composition for T=%f K'%T[i])
    _,comp,_ = HGS.eq(['H2','O2','H2O','H','O','OH'],[2,1,0,0,0,0],T[i],1)
    xH2[i]  = comp[0]/np.sum(comp)
    xO2[i]  = comp[1]/np.sum(comp)
    xH2O[i] = comp[2]/np.sum(comp)
    xH[i]   = comp[3]/np.sum(comp)
    xO[i]   = comp[4]/np.sum(comp)
    xOH[i]  = comp[5]/np.sum(comp)

# Print performance info
HGS.cr_info()

# plot the results and apply some basic formatting
plt.figure()
plt.plot(T,xH2O,'r',label='H2O',linewidth=2)
plt.plot(T,xH,  'b',label='H',  linewidth=2)
plt.plot(T,xO,  'g',label='O',  linewidth=2)
plt.plot(T,xOH, 'k',label='OH', linewidth=2)
plt.plot(T,xH2, 'c',label='H2', linewidth=2)
plt.plot(T,xO2, 'm',label='O2', linewidth=2)
plt.legend(loc='center left')
plt.xlabel('Temperature (K)')
plt.ylabel('Molar fraction')
plt.title(r'$H_2O$ dissociation as a function of temperature')
plt.grid()
plt.show()