#***********************************************************************************************************
# *HGSpy 
# *By Caleb Fuster, Manel Soria and Arnau Miró
# *ESEIAAT UPC      
#***********************************************************************************************************%
#
# Operating with mixtures: RP1 surrogate model and air
from __future__ import print_function

import numpy as np, HGSpy as HGS

HGS.set_options('warnings',False) # Deactivate warnings


# lets see what are RP1 and Air7921
print( HGS.rebuild('RP1',1,1) )
print( HGS.rebuild('air7921',1,1) )


# now RP1 can be used as a single species
species = ['RP1','O2','CO2','CO','H2O','OH','O','H']
n       = [1    , 1  , 0   , 0  , 0   , 0  , 0 , 0 ]
print( HGS.prop(species,n,1500,100) )

# RP1 and air combustion[
print( HGS.Tp(['RP1','air7921','CO2','CO','H2O','OH','O','H'],[1, 1, 0, 0, 0, 0, 0, 0 ],'T',300,40) ) # not O2 as it is already in Air7921 )

HGS.cr_info()