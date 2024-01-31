'''
***********************************************************************************************************
HGS CHEMICAL EQUATION SOLVER

HGS Solvers

By Caleb Fuster, Manel Soria and Arnau Miró
ESEIAAT UPC      
***********************************************************************************************************
'''
from __future__ import print_function, division

import numpy as np

from .cr import cr, cr_start, cr_stop


# Secant algorithm
#@cr('HGS.secant')
def hgs_secant(fun, n0, opt_sec, ilevel=0):
	"""
	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	species, n, Gmin = hgs_secant(fun, n0, opt_sec)

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*

	hgs_secant solves a function using a combination of the secant method and
	bisection method

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Inputs:
	-----------------------------------------------------------------------------
	fun --> Function
	n0 --> [mol] Initial mixture
	opt_sec --> Dictionary with the options for the secant method.
					"xmin": [K] Temperature minimum for the solver;
					"xmax" [K] Temperature maximum for the solver;
					"maxiter" Max iterations for the solver;
					"epsx" Diferential T where the solver reachs the solution;
					"epsy" Diferential S where the solver reachs the solution;
					"fchange" T difference where secant method is
							 changed by bisection method;
					"tipo" Select between: 'Frozen' for frozen flow
										  'Shifting' for shifting flow
					"info" Detailed info == 1; No info == 0.
					"dTp" Improve the velocity with the approximation of
						  parabola. +- dTp
			opt_sec = {"xmin": 300, "xmax": 4000, "maxiter": 200,
					   "epsx": 0.1, "epsy": 1, "tipo": "Shifting",
					   "fchange": 500, "info": 0, "dTp": 100}

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	Outputs:
	-----------------------------------------------------------------------------
	Tp --> [K] Final temperature
	n --> [mol] Final mixture
	flag --> Solver error detection:
				  1  Solver has reached the solution
				 -1  Solver failed. Maximum iterations
				 -2  Solver failed. Initial sign change not found

	*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*
	* Python HGS 1.0 from Matlab HGS 2.0
	* By Caleb Fuster, Manel Soria and Arnau Miró
	* ESEIAAT UPC
	"""
	cr_start('HGS.secant',ilevel)
	x1      = opt_sec['xmin']
	x2      = opt_sec['xmax']
	maxiter = opt_sec['maxiter']
	epsx    = opt_sec['epsx']
	epsy    = opt_sec['epsy']
	fchange = opt_sec['fchange']
	info    = opt_sec['info']
	dTp     = opt_sec['dTp']

	y1, n1 = fun(x1, n0)
	y2, n2 = fun(x2, n0)
	Tp, n  = [],[]

	if y1*y2 > 0:  # No sign change, sorry !
		flag = -2  # Initial sign change not found
		cr_stop('HGS.secant',ilevel)
		return Tp, n, flag

	if x2 - x1 > 1500:  # Try to fit to a parabola and solve the eq.
		x3 = (x2 + x1) / 2
		y3, _ = fun(x3, n0)
		a = (y1 - (y2 - y3) / (x2 - x3) * x1 - y3 +
			 x3 * (y2 - y3) / (x2 - x3)) / (x1 ** 2 + (x1 - x3) * (x3 ** 2 - x2 ** 2) / (x2 - x3) - x3 ** 2)
		b = (y2 - y3 + a * x3 ** 2 - a * x2 ** 2) / (x2 - x3)
		c = y3 - a * x3 ** 2 - b * x3

		xsol = [None, None]
		# Solution of the parabola
		xsol[0] = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
		xsol[1] = (-b - np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)

		if xsol[0] >= x1:
			x1p = xsol[0] - dTp
			x2p = xsol[0] + dTp
		else:
			x1p = xsol[1] - dTp
			x2p = xsol[1] + dTp

		# Don't allow the current range to extend the range specified
		if x1p < x1: x1p = x1
		if x2p > x2: x2p = x2

		y1p, n1 = fun(x1p, n0)
		y2p, n2 = fun(x2p, n0)

		if y1p*y2p < 0:  # Parabola method is okay
			x1 = x1p
			y1 = y1p
			x2 = x2p
			y2 = y2p

	flag  = -1  # We assume we are not solving it
	Bisec = False
	for ii in range(maxiter):

		if Bisec:  # If limits are close, switch to bisection algorithm
			xc = (x1 + x2) / 2
			Bisec = False
		else:
			xc = x1 - y1 * (x2 - x1) / (y2 - y1)  # Secant method

		if xc - x1 < x2 - xc:  # Next iteration is closer to x1
			n = n1
		else:
			n = n2

		yc, n = fun(xc, n)  # Compute next value

		if info:
			print(f'ii={ii} x1={x1:.2E} y1={y1:.2E} xc={xc:.2E} yc={yc:.2E} x2={x2:.2E} y2={y2:.2E} tolx={x2-x1:.2E}'
				  f' toly={yc:.2E}\n')

		if abs(yc) < epsy or (abs(xc - x1) < epsx and abs(x2 - xc) < epsx):  # Stop if it is solved
			flag = 1
			Tp = xc
			break

		if yc * y1 > 0:  # Change limits
			if x2 - xc < fchange or xc - x1 < fchange: Bisec = True
			y1 = yc
			x1 = xc
			n1 = n
		else:
			if x2 - xc < fchange or xc - x1 < fchange: Bisec = True
			y2 = yc
			x2 = xc
			n2 = n

	cr_stop('HGS.secant',ilevel)
	return Tp, n, flag