#!/usr/bin/env python
#
# HGSpy - HGS chemical equation solver in python
#
# By Caleb Fuster, Arnau Miro and Manel Soria
# UPC - ESEIAAT
from __future__ import print_function, division

import sys, os, numpy as np

from setuptools import setup, find_packages

with open('README.md') as f:
	readme = f.read()

# Main setup
setup(
	name             = 'HGSpy',
	version          = '1.1',
	author           = 'Caleb Fuster, Arnau Miro, Manel Soria',
	author_email     = 'calebfuji@gmail.com, arnau.miro@upc.edu, manel.soria@upc.edu',
	maintainer       = 'Arnau Miro',
	maintainer_email = 'arnau.miro@upc.edu',
	long_description = readme,
	url              = 'https://github.com/CalebFuster/HGSpy',
	packages         = find_packages(exclude=('Examples', 'doc')),
	install_requires = ['numpy','scipy','matplotlib']
)