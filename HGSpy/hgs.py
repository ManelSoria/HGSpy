'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS main class

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import numpy as np, pickle as pkl

from .            import HGSDATA
from .hgs_id      import hgs_id
from .hgs_mixture import hgs_add_mixture, hgs_subt_mixture, hgs_rebuild
from .hgs_print   import hgs_print_info
from .utils       import raiseError
from .definitions import R


class HGSData():
	'''
	'''
	def __init__(self, data={}):
		self._data = data

	def __len__(self):
		return len(self._data['name'])

	def __str__(self):
		return self._data.__str__()

	# Set and get functions
	def __getitem__(self,key):
		'''
		Recover the value of a variable given its key
		'''
		return self._data[key]

	def __setitem__(self,key,value):
		'''
		Set the value of a variable given its key
		'''
		self._data[key] = value

	# -- IDs --
	def id(self,species,raise_error=True):
		'''
		Run hgs_id
		'''
		return hgs_id(species,self._data,raise_error=raise_error)

	# -- Add, remove and rebuild --
	def add(self,name,species,percent):
		'''
		Run hgs_add_mixture
		'''
		return hgs_add_mixture(name,species,percent,hgs_data=self)

	def subt(self,name):
		'''
		Run hgs_subt_mixture
		'''
		return hgs_subt_mixture(name,self)

	def rebuild(self,species,n,T):
		'''
		Run hgs_rebuild
		'''
		return hgs_rebuild(species,n,T,self)

	def print_info(self,name):
		'''
		Run hgs_print_info
		'''
		hgs_print_info(name,self)

	# -- Properties --
	def coefs(self,ids,T,all=False):
		"""
		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		a = Coef_NASA(hgs_data, T, ids, species, **kwargs_coef)

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		Coef_NASA returns the NASA polynomials for a species at a certain temperature
		If they require both (low and high) values

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Inputs:
		-----------------------------------------------------------------------------
		T --> [K] Temperature
		id --> Sepecies id

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Outputs:
		-----------------------------------------------------------------------------
		a --> NASA polynomials

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		* Python HGS 1.0 from Matlab HGS 2.0
		* By Caleb Fuster, Manel Soria and Arnau Miró
		* ESEIAAT UPC
		"""
		if all:
			return self._data['lv'][ids], self._data['hv'][ids]

		lims = self._data['lim'][ids]
		if T < lims[0] or T > lims[2]:
			raiseError(f"hgs_single: Ups... Temperature {T} is not between the limits ({lims[0]:.2f}K-{lims[2]:.2f}K) for {self._data['name'][ids]}")

		return self._data['lv'][ids] if T <= lims[1] else self._data['hv'][ids]

	def cp(self,ids,T):
		"""
		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		cp = prop_cp(a, T)

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		prop_cp calculates the species Constant pressure coeficient using Burcat
		coeficients  and temperature

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Inputs:
		-----------------------------------------------------------------------------
		a --> Burcat coefficients
		T --> [K] Temperature

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Outputs:
		-----------------------------------------------------------------------------
		Cp --> [kJ/(mol*K)] Constant pressure coeficient

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		* Python HGS 1.0 from Matlab HGS 2.0
		* By Caleb Fuster, Manel Soria and Arnau Miró
		* ESEIAAT UPC
		"""
		a = self.coefs(ids,T)
		return R*(a[0] + np.sum([a[i]*T**i for i in range(1,5)])) # [kJ/(mol*K)]

	def cv(self,ids,T):
		"""
		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		cv = prop_cv(cp)

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		prop_cv calculates the species Constant volume coeficient using constant
		pressure coeficient

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Inputs:
		-----------------------------------------------------------------------------
		cp --> [kJ/(mol*K)] Constant pressure coeficient

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Outputs:
		-----------------------------------------------------------------------------
		cv --> [kJ/(mol*K)] Constant volume coefficient

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		* Python HGS 1.0 from Matlab HGS 2.0
		* By Caleb Fuster, Manel Soria and Arnau Miró
		* ESEIAAT UPC
		"""
		return self.cp(ids,T) - R # [kJ/(mol*K)]

	def h(self,ids,T):
		"""
		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		h = prop_h(a,T)

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		prop_h calculates the enthalpy of a species using his Burcat coeficients and
		temperature

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Inputs:
		-----------------------------------------------------------------------------
		a --> Burcat coefficients
		T --> [K] Temperature

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Outputs:
		-----------------------------------------------------------------------------
		h --> [kJ/mol] Enthalpy

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		* Python HGS 1.0 from Matlab HGS 2.0
		* By Caleb Fuster, Manel Soria and Arnau Miró
		* ESEIAAT UPC
		"""
		a = self.coefs(ids,T)
		return R*(a[5] + np.sum([a[i-1]*T**i/i for i in range(1,6)])) # [kJ/mol]

	def s(self,ids,T,P,Pref=1.):
		"""
		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		s = prop_s(a,T,P,state)

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		prop_s calculates the enthropy of a species using his Burcat coeficients
		and temperature

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Inputs:
		-----------------------------------------------------------------------------
		a --> Burcat coefficients
		T --> [K] Temperature
		P --> [bar] Pressure
		state -->  State of the species

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Outputs:
		-----------------------------------------------------------------------------
		s --> [kJ/(mol*K)] Entropy

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		* Python HGS 1.0 from Matlab HGS 2.0
		* By Caleb Fuster, Manel Soria and Arnau Miró
		* ESEIAAT UPC
		"""
		a = self.coefs(ids,T)
		s = R*(a[6] + a[0]*np.log(T) + np.sum([a[i]*T**i/i for i in range(1,5)]))  # [kJ/(mol*K)]
		if self._data['state'][ids] == "G" and not P == 0:
			s -= R*np.log(P/Pref) # [kJ/(mol*K)]
		return s

	def g(self,ids,T,P,Pref=1.):
		"""
		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		s = prop_s(s,h,T)

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

		prop_s calculates the species free Gibbs energy using enthalpy, enthropy
		and temperature

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Inputs:
		-----------------------------------------------------------------------------
		s --> [kJ/(mol*K)] Entropy
		h --> [kJ/mol] Enthalpy
		T --> [K] Temperature

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		Outputs:
		-----------------------------------------------------------------------------
		g --> [kJ/mol] Free Gibbs energy

		*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
		* Python HGS 1.0 from Matlab HGS 2.0
		* By Caleb Fuster, Manel Soria and Arnau Miró
		* ESEIAAT UPC
		"""
		return self.h(ids,T) - T*self.s(ids,T,P,Pref=Pref) # [kJ/mol]

	# -- Utilities --
	def append_dict(self,d):
		'''
		Append to the dataset
		'''
		for key in d.keys():
			self._data[key].append(d[key])

	def save(self,fname=HGSDATA):
		'''
		Save HGS data to file for exchange
		'''
		file = open(fname,'wb')
		pkl.dump(self._data,file)
		file.close()

	@classmethod
	def load(cls,fname=HGSDATA):
		'''
		Load HGS data
		'''
		file = open(fname,'rb')
		data = pkl.load(file)
		file.close()
		return cls(data=data)

	@classmethod
	def new(cls,name=[],nameback=[],state=[],lim=[],ena=[],nat=[],lv=[],
			hv=[],mm=[],comb=[],cspec=[],cper=[]):
		'''
		New empty class instance
		'''
		data = {
			'name'     : name,
			'nameback' : nameback,
			'state'    : state,
			'lim'      : lim,
			'ena'      : ena,
			'nat'      : nat,
			'lv'       : lv,
			'hv'       : hv,
			'mm'       : mm,
			'comb'     : comb,
			'cspec'    : cspec,
			'cper'     : cper,
		}
		return cls(data=data)