'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

Find the ID of species in the HGS database

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import numpy as np

from .cr    import cr, cr_stop
from .utils import raiseError


def find_hgs_id(name, hgs_data, combination, ns, raise_error):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	 ids = find_hgs_id(name, hgs_data, combination, ns)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	 ismember returns which index a string is found in a list of the database

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	 name --> String to find
	 hgs_data --> Database
	 combination --> Boolean to determinate if exist mixtures
	 ns --> Length of the database

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	 ids --> Index where the string was found

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	 * Python HGS 1.0 from Matlab HGS 2.0
	 * By Caleb Fuster, Manel Soria and Arnau Miró
	 * ESEIAAT UPC
	"""
	search  = np.array([name == s for s in hgs_data['name']])
	search1 = np.array([name == s for s in hgs_data['comb']] if combination else [0])

	if np.any(search) == 0 and np.any(search1) == 0:
		if raise_error:
			cr_stop('HGS.id',0)
			raiseError(f'hgs_id: {name} not found in the data base')
		else:
			return None
		
	ns = len(hgs_data['mm'])
	return np.where(search)[0][0] if np.any(search) else ns + np.where(search1)[0][0]


@cr('HGS.id')
def hgs_id(species,hgs_data,raise_error=True):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	ids = hgs_id(species,hgs_data)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_id finds the id of the species to improve the velocity of the code

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	species  --> String or numbers of species
	hgs_data --> HGS database

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	ids --> Species id

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	combination = not len(hgs_data['comb']) == 0

	if type(species) is list and all(type(s) is str for s in species):
		return [find_hgs_id(s,hgs_data,combination,len(hgs_data),raise_error) for s in species]

	if type(species) is str:
		return [find_hgs_id(species,hgs_data,combination,len(hgs_data),raise_error)]

	cr_stop('HGS.id',0)
	raiseError("uhh ? hgs_id wrong data type")