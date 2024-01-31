#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************
#
# CH4 combustion without dissociation
# 
# Inlet: 6.25 CH4 + 25 O2 + 25.79/21 N2 at Tr
# Outlet: 12.5 H2O + 12.5 O2 + 25.79/21 N2 at Tp
from __future__ import print_function

import numpy as np, HGSpy as HGS

from scipy.optimize import fsolve


Tr = 400 # K
p  = 1   # not rellevant 
 
# Specific enthalpy (h) of each species is obtained
hCH4 = lambda T : HGS.single('CH4','h',T,p)
hO2  = lambda T : HGS.single('O2' ,'h',T,p)
hN2  = lambda T : HGS.single('N2' ,'h',T,p)
hH2O = lambda T : HGS.single('H2O','h',T,p)
hCO2 = lambda T : HGS.single('CO2','h',T,p)
 
# Combustion equation to solve: sum(hreactives) - sum(hproducts) = 0
eq = lambda Tp : 6.25*hCH4(Tr) +   25*hO2(Tr) + 25*(79/21)*hN2(Tr) - \
                (12.5*hH2O(Tp) + 12.5*hO2(Tp) + 25*(79/21)*hN2(Tp) + 6.25*hCO2(Tp))

Tp = fsolve(eq,1000.,xtol=1e-4,maxfev=4000)
print('Adiabatic flame temperature: %.3f K'%Tp)

# Print performance info
HGS.cr_info()