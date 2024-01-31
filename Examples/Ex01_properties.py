#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************
#
#  Properties of single elements and mixtures
from __future__ import print_function

import numpy as np, HGSpy as HGS


Ru = HGS.R # kJ/molK

# Compute the propeties of a mixture of 1 mol of N2, 2 mol of CH4 and 3 mol
# of C3H8 at 3 bar and 400K 
Mm,Cp,Cv,H,S,G,Rg,gamma,a = HGS.prop(['N2','CH4','C3H8'],[1,2,3],400,3)
print('propeties of a mixture of 1 mol of N2 , 2 mol of CH4 and 3 mol of C3H8 at 3 bar and 400K')
print('Mm = %f, Cp = %f, Cv = %f'%(Mm,Cp,Cv))
print()

# Using N2, verify that deltaH is aprox. equal to Cp*deltaT for small deltaT
CpN2,hN2_1 = HGS.prop('N2',1,300,10,'Cp','H')
hN2_2      = HGS.prop('N2',1,301,10,'H')
print('Single species  N2')
print('DeltaH = <%.4f> , Cp = <%.4f>'%((hN2_2-hN2_1),CpN2))
print()

Cpmix,hmix_1 = HGS.prop(['N2','O2'],[1,1],300,10,'Cp','H');
hmix_2       = HGS.prop(['N2','O2'],[1,1],301,10,'H');
print('Mixture of N2 and O2')
print('DeltaH = <%f> , Cp = <%f>'%((hmix_2-hmix_1),Cpmix))
print()

# Using N2, verify that S2-S1 aprox. = cp.ln(T2/T1)-Ru*ln(P2/P1)
S1          = HGS.prop('N2',1,400,20,'S')
CpN2,RN2,S2 = HGS.prop('N2',1,500,10,'Cp','Rg','S')
print('Single species  N2')
print('S2-S1 = <%f> , cp*ln(T2/T1)-Ru*ln(P2/P1) = <%f>'%(S2-S1,CpN2*np.log(500/400)-Ru*np.log(10/20)))
print()

# Using air (80% N2, 20% O2) 
# verify @300K deltaH aprox. equal to Cp*deltaT
CpA,hA_1 = HGS.prop(['N2','O2'],[0.8,0.2],300,10,'Cp','H')
hA_2     = HGS.prop(['N2','O2'],[0.8,0.2],310,10,'H')
print('Mixture of N2 and O2')
print('dH/10 = <%f> , Cp = <%f>'%((hA_2-hA_1)/10,CpA))
print()

# @300K and 1 bar verify:
T = 300 # K
P = 1   # bar
m = 10  # kg (random)

MM,Cp,Cv,H,S,G,Rg,gamma,a = HGS.prop(['N2','O2'],[0.8,0.2],T,P)

# sound speed  (J/kgK)^(1/2)=m/s
print('Mixture of N2 and O2')
print('a = <%f> , sqrt(gamma*Rg*1000*T) = <%f> '%(a,np.sqrt(gamma*Rg*1000*T)))

# Rg=Ru/MM  kJ/kg K
# The code is printing here the gas coefficient R
print('1000*Ru/MM = <%f> , Rg = <%f>'%(1000*Ru/MM ,Rg))

# Cp/Cv=gamma
print('Cp/Cv-gamma = <%f>'%(Cp/Cv-gamma))

# g=h-T*s
v = 1000*Rg*T/(P*1e5) # m^3/kg
V = m*v
n = m*1000/MM # mol

h = H/m # kJ/kg
s = S/m # kJ/kgK
g = G/m # kJ/kg
u = h-P*1e5*v/1000
print('h-T*s = <%f>,g = <%f>'%(h-T*s,g))
print()


# verify molar cp of the mixture
CpO2 = HGS.single('O2','cp',T,1) # kJ/molK
MMO2 = HGS.single('O2','Mm',T,1) 
CpN2 = HGS.single('N2','cp',T,1) # kJ/molK
MMN2 = HGS.single('N2','Mm',T,1) 
hN2  = HGS.single('N2','h', T,1) # kJ/mol
print('Verifying independently with HGSsingle (0.8*CpN2+0.2*CpO2-Cp) = <%f>'%(0.8*CpN2+0.2*CpO2-Cp))

CpN2kg = 1000*CpN2/MMN2 # kJ/molK -> kJ/kgK
hN2kg  = 1000*hN2/MMN2  # kJ/mol -> kJ/kgK 
print('N2: Cp=%f kJ/kgK hN2 =%f kJ/kg'%(CpN2kg,hN2kg))

# Note that the program returns in units of kJ/molK. Here are converted to
# kJ/kgK
Cpkg = Cp/(MM/1000)
print('cp=%f kJ/kgK'%Cpkg)
print()

# hgssingle computes H, G and S properties for a single element
print('hgssingle computes H, G and S properties for a single element')
print( HGS.single('O2','h',500,10) )

# Print performance info
HGS.cr_info()
