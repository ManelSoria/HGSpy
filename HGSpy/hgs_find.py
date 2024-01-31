'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS Mixture main functions

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

from .cr          import cr


@cr('HGS.find')
def hgs_find(name,complete,hgs_data):
    '''
    **************************************************************************
    
     HGSfind(name,complete)
    
    **************************************************************************
    
     HGSfind finds the species and the mixtures that contain the string name.
    
    **************************************************************************
     Inputs:
    --------------------------------------------------------------------------
     name --> String to be found
     complete --> Prints the complete name if complete = 1
    
     Outputs:
    --------------------------------------------------------------------------
     Command Window Print
    
    **************************************************************************
    * Python HGS 1.0 from Matlab HGS 2.1
    * By Caleb Fuster, Manel Soria and Arnau Miró
    * ESEIAAT UPC 
    '''
    name_base = hgs_data['name'] if not complete else hgs_data['nameback']

    # Search species
    cap = True # in case there is no name with the string
    print('Species that contain %s'%name)
    for ii,n in enumerate(name_base):
        if name in n:
            cap = False
            print('<%d>  %s'%(ii,n))
    if cap: # Print None to be a fancy list
        print('None')

    # Search combinations
    cap = True
    print('Mixtures that contain %s'%name)
    for ii,n in enumerate(hgs_data['comb']):
        if name in n:
           cap = False
           print('<%d>  %s'%(len(hgs_data)+ii,n))
    if cap: # Print None to be a fancy list
        print('None')