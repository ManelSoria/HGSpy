#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************
#
# Turbopump example (liquid H2O2 adiabatic decomposition)
#
# 80% H2O2 mass fraction
# vaporization heat of H2O2 @ 80%: 420 cal/g = 420*4.18 kJ/kg
# 
# Inlet:  molH2O2.H2O2 + molH2O.H2O
# Outlet: mollH2O2. (H2O + (1/2)O2 ) + molH2O a Tp
from __future__ import print_function

import numpy as np, HGSpy as HGS

from scipy.optimize import fsolve
 

R  = HGS.R # kJ/molK
Tr = 300   # K
pp = 8     # bar (products' pressure)

# We use 1Kg as a base for our calculations
molH2O2 = 0.8*(1/34e-3) # molH2O2 for kg of reactives
molH2O  = 0.2*(1/18e-3) # molH2O  for kg of reactives

# Computation of enthalpy, pressure is not rellevant
hH2O2 = lambda T: HGS.single('H2O2','h',T,1)
hO2   = lambda T: HGS.single('O2',  'h',T,1)
hH2O  = lambda T: HGS.single('H2O', 'h',T,1)

# Equation to be solved
eq = lambda Tp: molH2O2*hH2O2(Tr) + molH2O*hH2O(Tr) - 420*4.18 \
              -(molH2O2*(hH2O(Tp) + 0.5*hO2(Tp))    + molH2O*hH2O(Tp))

Tp = fsolve(eq,3000.,xtol=1e-4,maxfev=4000)
print('Tp = %.2f K'%Tp)

# Computation of the turbopump power
# first, gamma coeficient is needed
Cp,MM,gamma = HGS.prop(['H2O','O2'],[molH2O2+molH2O,molH2O2*0.5],Tp,pp,'Cp','Mm','gamma')

p2  = 1 # pressure at turbine outlet
nus = 0.85

T2s    = Tp*(p2/pp)**( (gamma-1)/gamma )
deltaT = nus*(Tp - T2s)
T2     = Tp - deltaT

# Cp is in kJ/molK, and needs to be converted to kJ/kgK:
# MM in g/mol
Cp = Cp/(MM*1e-3)
w  = Cp*deltaT # kJ/kg of propulsant
print(Cp, w)

# Print performance info
HGS.cr_info()