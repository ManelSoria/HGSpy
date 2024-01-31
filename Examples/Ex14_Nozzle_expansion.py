#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************%
#
# Nozzle expansion process.
# LH2-LOX reaction
#
# Exemplification of the differents variables across of a nozzle expansion
# Changing the for of the P vector by a solver is what HGSnozzle does.
from __future__ import print_function

import numpy as np, HGSpy as HGS
import matplotlib.pyplot as plt


HGS.set_options('warnings',False) # Deactivate warnings

# The inlet conditions to combustion chamber of the Vinci engine are:
species = [    'H2',    'O2','H2O','OH','O','H']
n0      = [2875.496,1052.566,    0,   0,  0,  0] # mol/s (obtained from kg/s)
P0      = 62
Pa      = 0        # Vacuum
Hin     = -19.0609 # kJ (liquid inlet, obtained with INIST)

Tp,npp,_,_ = HGS.Tp(species,n0,'H',Hin,P0)


## Nozzle expansion
# We generate a vector with the pressure points to be obtained
# From 10 bar below the chamber pressure to 0.01 in 2 diferent steps
# 10 bar is arbitrary, we just want to make sure that we begin before the
# throat but far from the inlet, where the velocity would be 0 and the area
# infinite.
P = np.arange(58,1,-2).tolist() + np.linspace(1,0.01,20).tolist()

species,n,T,v,M,A,F,Isp = HGS.nozzle(species,npp,Tp,P0,P,Pa)


# Plots
plt.figure()
plt.plot(P,T,'-b',linewidth=2)
ax = plt.gca()
ax.invert_xaxis()
plt.title('T versus P')
plt.xlabel('P [bar]')
plt.ylabel('T [K]')

plt.figure()
plt.plot(P,n[2,:]/np.sum(n[2,:]),linewidth=2)
ax = plt.gca()
ax.invert_xaxis()
plt.title('H2O fraction versus P')
plt.xlabel('P [bar]')
plt.ylabel('H2O Molar fraction []')

plt.figure()
plt.plot(P,v,linewidth=2)
plt.title('v vs P')
ax = plt.gca()
ax.invert_xaxis()
plt.xlabel('P [bar]')
plt.ylabel('v [m/s]')

plt.figure()
plt.plot(P,M,linewidth=2)
plt.plot(P,np.ones_like(P),'r')
plt.title('Mach vs P')
ax = plt.gca()
ax.invert_xaxis()
plt.xlabel('P [bar]')
plt.ylabel('Mach []')
# We can find the exact throat pressure with:
_,_,_,_,Pt,flag = HGS.isentropic(species,npp,Tp,P0,'M',1)
plt.plot([Pt,Pt],[0,6],'r')

plt.figure()
plt.plot(P,F,linewidth=2)
plt.title('Thrust vs P')
ax = plt.gca()
ax.invert_xaxis()
plt.xlabel('P [bar]')
plt.ylabel('F [N]')

plt.figure()
D = 1000*2*np.sqrt(A/np.pi)
plt.semilogy(P,D,linewidth=2)
plt.semilogy([Pt,Pt],[np.min(D)/2,np.max(D)],'r')
plt.title('Diameter vs P')
plt.xlabel('P [bar]')
plt.ylabel('Diameter [mm]')

Dt = np.interp(P0-Pt,P0-np.array(P),D) # P0-P to have x values in ascending order
De = D[-1]

print('Throat diameter %.2f mm\nExit diameter %.2f mm'%(Dt,De))
print('Nozzle expansion ratio = %.2f'%((De/Dt)**2))
print('Thrust = %.2f kN'%(F[-1]/1000))

HGS.cr_info()
plt.show()