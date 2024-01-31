#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************%
#
# Example by Caleb Fuster
# The MAN problem
# Wet air reaction with hydrocarbons. Using Relative humidity.
# At ambient pressure and temperature.
from __future__ import print_function

import numpy as np, HGSpy as HGS
import matplotlib.pyplot as plt


HGS.set_options('warnings',False) # Deactivate warnings

RHv      = np.arange(0,100,2) # [%] Relative humidity
T1       = 298.15             # [K] Air temperature
T2       = 350                # [K] Fuel temperature
Kgwetair = 1                  # [kg] Wet air mass flow
P        = 1                  # [bar] Atmospheric pressure
V        = 0.01               # [m^3] Volume
R        = 8.3144621          # [J / (mol*K)]

nwetair  = P*10**5*V/(R*T1)   # [mol]

Tp = np.zeros_like(RHv)

opt_sec = HGS.options
opt_sec['xmin'] = 500
opt_sec['xmax'] = 2500
opt_sec['epsy'] = 1
opt_sec['epsx'] = 2

opt_eq  = {
    'method':'SLSQP',
    'tol':1e-6,
    'options':{},
}

# Good approx up to 40ºC. Be carefull Tc (ºC)
# http://hyperphysics.phy-astr.gsu.edu/hbase/Kinetic/relhum.html#c4
VDsat    = lambda Tc : 5.018 + 0.32321*Tc + 8.1847e-3*Tc**2 +3.1243e-4*Tc**3
VDsatAmb = VDsat(T1-273.15)

# Dry air composition in mols
dryaircomp = [0.79,0.21];

# Molar mass
Mmdryair = HGS.prop(['N2','O2'],dryaircomp,1,1,'Mm')[0]
MmH2O    = HGS.prop(['H2O'],1,1,1,'Mm')[0]

# Wet air
wet = ['N2','O2','H2O'];

# Fuels
fuel  = ['CH4','C2H6','C3H8']
nfuel = [1.0,0.8,0.5]  # [mols]

# Products
prod  = ['CO2','CO']
nprod = [0,0]

n    = np.zeros((len(RHv),len(wet+fuel+prod)))
nper = np.zeros_like(n)
for ii,RH in enumerate(RHv):
    print('Running for RH=%.2f%%'%(RH))

    # RH = VD / VDsat *100;
    VD = RH*VDsatAmb/100
    
    # VD = KgH2O / KgH2 =  KgH2O / 2
    KgH2O = VD*2
    nH2O  = KgH2O/MmH2O
            
    ndryair  = nwetair - nH2O
    
    # Wetair composition in mols
    nwetaircomp = [dryaircomp[0]*ndryair,dryaircomp[1]*ndryair,nH2O]/(ndryair + nH2O)

    # Wet air
    nwetair = Kgwetair*1000/HGS.prop(wet,nwetaircomp,1,1,'Mm')[0] # [mols] 
    nwet    =  (nwetair*nwetaircomp).tolist() # [mols]

    # Species
    # HGS functions requires all the species that you want to be calculated.
    # If you want to put dissociation, it is possible. Like in this example H2
    # and H or O2 and O.
    species = wet  + fuel  + prod
    n0      = nwet + nfuel + nprod

    # Assuming that Wet air and Fuels have not the same inlet T
    Hwet   = HGS.prop(wet,nwet,T1,1,'H')[0]
    Hfuels = HGS.prop(fuel,nfuel,T2,1,'H')[0]
    Hin    = Hwet + Hfuels;

    # HGStp
    Tp[ii],n[ii,:],species,_ = HGS.Tp(species,n0,'H',Hin,1,opt_eq=opt_eq,opt_sec=opt_sec)
    nper[ii,:] = n[ii,:]/np.sum(n[ii,:])

plt.figure()
plt.plot(RHv,Tp)
plt.ylabel('Tp [K]')
plt.xlabel('Humidity [%]')

plt.figure()
for ii,s in enumerate(species):
    plt.plot(RHv,nper[:,ii]*100,label=s)
plt.xlabel('Humidity [%]')
plt.ylabel('Mols [%]')
plt.axis([0,100,0,25])
plt.legend()

plt.figure()
for ii,s in enumerate(species):
    plt.plot(RHv,n[:,ii],label=s)
plt.xlabel('Humidity [%]')
plt.ylabel('Mols [mols]')
plt.axis([0,100,0,5])
plt.legend()

HGS.cr_info()
plt.show()