#######################################
# spin-up Hyun and Park JFM, 1992;
# A compact FDM on staggered grid for Navier-Stokes flows, Zhang et. al., Int. J. Numer. Meth. Fluids. 2006;
# Author: Pratik Aghor
import numpy as np
from numpy.linalg import *

#######################################
#######################################
# will be using QR factorization to solve Ax = b many times
# so it makes sense to have a function for that.
def qrsolve(A, b):
	"""
	solve Ax = b using QR
	"""
	Q, R = qr(A)

	# solve Q*temp = b
	temp = np.matmul(np.transpose(Q), b)

	# solve R*x = temp
	x = solve(R, temp)

	return x
#######################################
#######################################
