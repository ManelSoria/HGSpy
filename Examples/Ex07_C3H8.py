#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************%
# Adiabatic C3H8 reaction using hgsTp
#
# Inlet: H2, O2 
# Outlet: H2O + (1/2)O2 at Tp
# TBD: rewrite the example so that the variable is:
# (1) the OF ratio 
# (2) the O2 excess
from __future__ import print_function

import numpy as np, HGSpy as HGS
import matplotlib.pyplot as plt


HGS.set_options('warnings',False) # Deactivate warnings

species =[
    'C3H8',
    'CO2',
    'CO',
    'O2',
    'O',
    'H2',
    'H',
    'OH',
    'H2O'
]

molprop = np.linspace(1,4,15)
Tp      = np.zeros_like(molprop)
for i in range(len(molprop)):
    print('Solving for i=%d / %d ... '%(i+1,len(molprop)))
    nr = [molprop[i],0,0,5,0,0,0,0,0] # mol
    Tp[i],_,np,_ = HGS.Tp(species,nr,'T',298,1)
HGS.cr_info() # Print performance info

plt.figure()
plt.plot(molprop,Tp,linewidth=2)
plt.grid();
plt.xlabel('n (mol C3H8)')
plt.ylabel('T (K)')
plt.title('Flame temperature (with dissociation) of n mol C3H8 + 5 mol O2')
plt.show()