'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS Mixture main functions

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

from .utils import raiseError


def hgs_print_info(name,hgs_data):
    """
    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    hgs_print_info(name, **kwargs)

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

    hgs_print_info finds the species and the mixtures that contains the (name) in
    his name.

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Inputs:
    -----------------------------------------------------------------------------
    name --> String of species
    **kwargs --> complete= to enable the full name return

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    Outputs:
    -----------------------------------------------------------------------------

    *+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
    * Python HGS 1.0 from Matlab HGS 2.0
    * By Caleb Fuster, Manel Soria and Arnau Miró
    * ESEIAAT UPC
    """
    ids = hgs_data.id(name)

    if not len(ids) == 1:
        raiseError("Please call HGSprintInfo with just a component")
    ids = ids[0]

    ns = len(hgs_data)
    if ns > ids:
        print(f"Species = <{hgs_data['name'][ids]}>   code = {ids}", end="\n")
        ena = hgs_data["ena"][ids]  # element names
        nat = hgs_data["nat"][ids]
        ne  = len(ena)  # number of elements
        print("- Composition: ", end=" ")
        for jj in range(ne):
            print(f"{ena[jj]}-{nat[jj]:.0f}  ", end=" ")
        print(f"\n- Mm = {hgs_data['mm'][ids]:.4f} ", end="\n")
        print(' -----------------\n')
    else:
        ids -= ns
        print(f'Combination = <{hgs_data["comb"][ids]}>   code = {ids + ns}', end="\n")
        cen = hgs_data["cspec"][ids] # combination element names
        cna = hgs_data["cper"][ids]  # combination element atoms
        ne  = len(cna)               # number of  elements
        print('- Composition  ', end=" ")
        for jj in range(ne):
            print(f'{cen[jj]}: { cna[jj]:.2f}%.  ', end=" ")
        print('\n -----------------\n')