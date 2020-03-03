import numpy as np  # Import numpy
import math
from numpy.linalg import *
#######################################
'''
Code to integrate 1d heat equation (u_t = nu*u_xx) with Dirichlet boundaries
using finite difference method. Crank-Nicholson CN was used for time marching.
Author: Pratik Aghor

This file contains required functions - to be called from main

# I => global grid, including ghost nodes
'''
#######################################
#######################################
def set_periodic_BCs(fI):
    n = len(fI); # I => global grid, including ghost nodes
    fI[0] = fI[n-2];
    fI[n-1] = fI[1];
    return fI
#######################################
#######################################
def getfI_1x(fI, dx):
    """
    first derivative of fI, periodic BC
    """
    n = len(fI)
    fI_1x = np.zeros(n)
    for j in range (1, n-1):
        fI_1x[j] = (fI[j+1]-fI[j-1])/(2.0*dx);

    set_periodic_BCs(fI);

    return fI_1x
#######################################
#######################################
def get_Burger_rhs(uI, dx):
    """
    RHS of Burger's equation = -uu_x
    """
    n = len(uI);
    rhs = np.zeros(n);
    set_periodic_BCs(uI);

    uI_1x = getfI_1x(uI, dx)
    # set_periodic_BCs(uI); already done inside getfI_1x

    rhs = -np.multiply(uI, uI_1x);
    return rhs
#######################################
#######################################
