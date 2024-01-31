'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

Data downloading routines

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import os, urllib.request, numpy as np

from ..            import RAWDATA
from ..hgs         import HGSData
from ..cr          import cr
from ..utils       import raiseError
from ..definitions import Mendeley


URLDATA = 'http://garfield.chem.elte.hu/Burcat/THERM.DAT'


def fix_str_float(s):
	return s.replace('D','E').replace('E ','E+').replace(' ','')


def parse_chemkin_elem(e):
	'''
	Parse a single element of the database
	'''
	out = {}
	# Parse name up to the first whitespace
	i1 = e.find(' ')
	n = e[:i1].strip()
	if n.lower() == 'air':   return None # Skip air as no composition is given
	if 'mgcl2' in n.lower(): return None # Skip as it is horrible to parse
	out['name'] = n
	e = e[18+6:]

	# Then parse composition, with a maximum
	# of 4 elements
	ename, enum = [], []
	for jj in range(4):
		# for that find the jth point
		# and parse the number
		ee = e[:5]
		if len(ee.strip()) == 0: e = e[5:]; continue
		n  = int(float(ee[2:].strip()))
		# If the number is 0 we do not
		# have more elements
		if n == 0: e = e[5:]; continue
		# Now we do have an element
		enum.append(n)
		# An element name has a maximum of two letters
		n = ee[:2].strip()
		# If len(n) == 2 make the second letter lower case
		if len(n) == 2: n = n[0] + n[1].lower()
		ename.append(n)
		e = e[5:]
	out['ena'] = ename
	out['nat'] = enum

	# Parse state
	out['state'] = e[:4].strip()
	e = e[4:]

	# Parse limits
	out['lim'] = [float(s.strip().replace(',','.')) for s in e[:25].split()]
	l = out['lim'][1]
	out['lim'][1] = out['lim'][2]
	out['lim'][2] = l
	e = e[e.find('\n')+1:]

	# Parse coefficients
	coefs = []
	for ii in range(3): # Read 3 rows of coefficients
		l = e[:e.find('\n')]
		for jj in range(5): # with 5 coefficients each
			coefs.append( float(fix_str_float(l[:15])) if not 'N/A' in l[:15] else np.nan )
			l  = l[15:]
		e = e[e.find('\n')+1:]
	coefs = np.array(coefs)
	out['hv'] = coefs[:7]
	out['lv'] = coefs[7:14]

	return out


@cr('HGS.download')
def download(download=True,info=False):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	 hgs_data_download(**kwargs)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	 downloads HGS database

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	 **kwargs --> download= To download the database from internet
				  info= To show the evolution of the database building

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	 * Python HGS 1.0 from Matlab HGS 2.0
	 * By Caleb Fuster, Manel Soria and Arnau Miró
	 * ESEIAAT UPC
	"""
	# ----- Downloading database ------ #
	if download:
		print("Downloading HGS database - \n")
		# URL
		urllib.request.urlretrieve(URLDATA,RAWDATA)
		print("HGS database downloaded - \n")
	if not os.path.isfile(RAWDATA): raiseError('Database not expected in the right directory <%s>!'%RAWDATA)

	try:
		hgs_old = HGSData.load(fname=HGSDATA)
		hgs_new = HGSData.new(comb=hgs_old['comb'],cspec=hgs_old['cspec'],cper=hgs_old['cper'])
	except:
		hgs_new = HGSData.new()

	# ----- Database processing ------ #
	# Open file read and text mode
	f = open(RAWDATA, 'r')
	lines = f.readlines()
	f.close()

	# Find index of THERMO ALL to start parsing
	# the database
	if info: print('Skipping header...')
	iline = -1
	for ii,l in enumerate(lines):
		if 'THERMO ALL' == l.strip():
			iline = ii
			break

	# Group the items in 4 by 4 so that
	# it corresponds to a full entry of the
	# database
	if info: print('Grouping table items...')
	dbase = []
	for ii in range((len(lines)-iline-1)//4):
		s, istart = '', iline + 2 + 4*ii
		s += lines[istart]
		for jj in range(3):
			s += lines[istart+jj+1]
		dbase.append(s)

	# Parse database items one by one
	if info: print('Parsing database...',end=' ')
	for d in dbase:
		out = parse_chemkin_elem(d)
		if out is not None: hgs_new.append_dict(out)
	if info: print('Done!')
	
	# Perform some name modifications and cleaning up
	# of the database
	if info: print('Name modifications...',end=' ')
	for ii in range(len(hgs_new)):
		# Fix that some names do not contain the element
		# by its lower case value
		for s in hgs_new['ena'][ii]:
			if s.upper() in hgs_new['name'][ii]: 
				hgs_new['name'][ii] = hgs_new['name'][ii].replace(s.upper(),s)
		
		# - Liquid Change - #
		for s in ['(L)','(liq)']:
			if s in hgs_new['name'][ii]:
				hgs_new['name'][ii] = hgs_new['name'][ii].replace(s,'(l)')

		# - Solid Change - #
		if '(S)' in hgs_new['name'][ii]:
			hgs_new['name'][ii] = hgs_new['name'][ii].replace('(S)','(s)')

		# Solid change  cr -> (cr)
		if ' cr' in hgs_new['name'][ii]:
			hgs_new['name'][ii] = hgs_new['name'][ii].replace(' cr','(cr)')
		if '(cr)A' in hgs_new['name'][ii]:
			hgs_new['name'][ii] = hgs_new['name'][ii].replace('(cr)A','(cr.A)')
		if '(cr)B' in hgs_new['name'][ii]:
			hgs_new['name'][ii] = hgs_new['name'][ii].replace('(cr)B','(cr.B)')
		if '(cr)C' in hgs_new['name'][ii]:
			hgs_new['name'][ii] = hgs_new['name'][ii].replace('(cr)C','(cr.C)')
		if '(cr)I' in hgs_new['name'][ii]:
			hgs_new['name'][ii] = hgs_new['name'][ii].replace('(cr)I','(cr.I)')
		if '(cr)II' in hgs_new['name'][ii]:
			hgs_new['name'][ii] = hgs_new['name'][ii].replace('(cr)II','(cr.II)') 
		if '(III)cr' in hgs_new['name'][ii]:
			hgs_new['name'][ii] = hgs_new['name'][ii].replace('(III)cr','(cr.III)')       

		# - Gas - #
		for s in ['(G)','(gas)',' gas',' g ']:
			if s in hgs_new['name'][ii]:
				hgs_new['name'][ii] = hgs_new['name'][ii].replace(s,'(g)')
		
		# - Ref element - #
		for s in ['REF','REF-ELEMENT']:
			if s in hgs_new['name'][ii]:
				hgs_new['name'][ii] = hgs_new['name'][ii].replace(s,'')

	if info: print('Done!')

	# Molar Mass calculation
	if info: print('Molar Mass calculation ...', end=' ')
	for ii in range(len(hgs_new)):
		mm = 0.
		for n,e in zip(hgs_new['nat'][ii],hgs_new['ena'][ii]):
			if e in Mendeley.keys(): mm += n*Mendeley[e]
		hgs_new["mm"].append(mm)
	if info: print('Done!')

	return hgs_new