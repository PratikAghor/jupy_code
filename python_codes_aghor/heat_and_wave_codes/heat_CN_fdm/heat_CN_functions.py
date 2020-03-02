import numpy as np  # Import numpy
import math
from numpy.linalg import *
#######################################
'''
Code to integrate 1d heat equation (u_t = nu*u_xx) with Dirichlet boundaries
using finite difference method. Crank-Nicholson CN was used for time marching.
Author: Pratik Aghor

This file contains required functions - to be called from main
'''
#######################################
def build_A(Nx, alpha):
    """
    function to build the LHS matrix in the CN formulation
    """
    # declare A
    A = np.zeros((Nx, Nx));

    # start building A, first and last rows represent Dirichlet BC's

    A[0, 0] = 1.0;
    for i in range(1, Nx-1):
        A[i, i-1] = -alpha;
        A[i, i]   = (1.0 + 2.0*alpha);
        A[i, i+1] = -alpha;
    A[Nx-1, Nx-1] = 1.0;

    return A
#######################################
def build_B(Nx, alpha):
    """
    function to build the RHS matrix in the CN formulation
    """
    # declare B
    B = np.zeros((Nx, Nx));

    # start building B, first and last rows represent Dirichlet BC's

    B[0, 0] = 1.0;
    for i in range(1, Nx-1):
        B[i, i-1] = alpha;
        B[i, i]   = (1.0 - 2.0*alpha);
        B[i, i+1] = alpha;
    B[Nx-1, Nx-1] = 1.0;

    return B
#######################################
def set_Dirichlet_BCs(u, ul, ur):
    u[0] = ul;
    u[-1] = ur;
    return u
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
