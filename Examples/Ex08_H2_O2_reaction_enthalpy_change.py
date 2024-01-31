#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************%
# Adiabatic H2 / O2 reaction with dissociation
# (as an example of how HPStp looks for the temperature)
#
# Inlet: H2, O2 
# Outlet: H2O + (1/2)O2 at Tp
from __future__ import print_function

import numpy as np, HGSpy as HGS
import matplotlib.pyplot as plt


species = ['H2','O2','H2O','H','O','OH']
nr      = [   2,   1,    0,  0,  0,   0] # mol
Tr      = 350 # K 
P       = 66  # bar
 
Hin     = HGS.prop(species,nr,Tr,P,'H')
 
T    = np.linspace(300,5000,20)
Hout = np.zeros_like(T)
for i in range(len(T)): 
    _,comp,_ = HGS.eq(species,nr,T[i],P)
    Hout[i]  = HGS.prop(species,comp,T[i],P,'H')[0]
 
plt.plot(T,Hout,'-r',label='Hout',linewidth=2)
plt.plot(T,Hin*np.ones_like(T),'-b',label='Hin',linewidth=2)


Tp,npp,_,_ = HGS.Tp(species,nr,'T',Tr,P)
print(npp/np.sum(npp))
flameHout = HGS.prop(species,npp,Tp,P,'H')
plt.plot(Tp,flameHout,'ok',label='HGStp reaction temperature',markersize=14,markerfacecolor='k')
plt.xlabel('T')
plt.ylabel('H')
plt.legend()
plt.title('H2 O2 combustion with dissociation')

HGS.cr_info()
plt.show()