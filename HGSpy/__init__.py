'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
import os

# Paths to important files
DATAPATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),'data')
HGSDATA  = os.path.join(DATAPATH,'data.hgs')
RAWDATA  = os.path.join(DATAPATH,'DATA_7_coef.txt')

from .            import data
from .hgs         import HGSData
from .definitions import R
from .cr          import cr_start, cr_stop, cr_reset, cr_info
from .utils       import set_options, get_options

# HGS functions
from .hgs_id         import hgs_id
from .hgs_find       import hgs_find
from .hgs_mixture    import hgs_add_mixture, hgs_subt_mixture, hgs_rebuild
from .hgs_print      import hgs_print_info
from .hgs_prop       import hgs_prop as prop, hgs_single as single
from .hgs_eq         import hgs_eq as eq
from .hgs_Tp         import hgs_Tp as Tp
from .hgs_isentropic import hgs_isentropic as isentropic
from .hgs_nozzle     import hgs_nozzle as nozzle
from .hgs_solver     import options

# Some predefined functions
id           = lambda species,hgs_data=HGSData.load(),raise_error=True : hgs_id(species,hgs_data,raise_error)
add_mixture  = lambda name,species,percent,hgs_data=HGSData.load()     : hgs_add_mixture(name,species,percent,hgs_data)
subt_mixture = lambda name,hgs_data=HGSData.load()                     : hgs_subt_mixture(name,hgs_data)
rebuild      = lambda species,n,T,hgs_data=HGSData.load()              : hgs_rebuild(species,n,T,hgs_data)
print_info   = lambda name,hgs_data=HGSData.load()                     : hgs_print_info(name,hgs_data)
find         = lambda name,complete=False,hgs_data=HGSData.load()      : hgs_find(name,complete,hgs_data)

del os, hgs, definitions, cr, utils
del hgs_id, hgs_prop, hgs_solver, hgs_mixture, hgs_print
del hgs_eq, hgs_Tp, hgs_isentropic, hgs_nozzle