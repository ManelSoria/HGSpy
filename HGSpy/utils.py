'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

Utilities functions, not callable within HGS frame

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import sys, numpy as np

OPTS = {'warnings':True,'errors':True}


def set_options(key,value):
	'''
	Set general HGSpy options
	'''
	OPTS[key.lower()] = value

def get_options(key):
	'''
	Get general HGSpy options
	'''
	return OPTS[key.lower()]


def raiseError(errmsg):
	'''
	Raise a controlled error and abort execution on
	all processes.
	'''
	if OPTS['errors']:
		print('Error: %s' % (errmsg),file=sys.stderr,flush=True)
		sys.exit(1)


def raiseWarning(warnmsg):
	'''
	Raise a controlled warning but don't abort execution on
	all processes.
	'''
	if OPTS['warnings']: print('Warning: %s' % (warnmsg),file=sys.stderr,flush=True)


def truncate(value,precision):
	'''
	Truncate array by a certain precision
	'''
	fact  = 10**precision
	return np.round(value*fact)/fact
