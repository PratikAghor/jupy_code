#######################################
# Computation of Strongly Compressible Rotating Flows, Harada, J. Comp. Phys., 1980;
# A compact FDM on staggered grid for Navier-Stokes flows, Zhang et. al., Int. J. Numer. Meth. Fluids. 2006;
# Author: Pratik Aghor
#######################################
import numpy as np
from numpy.linalg import *

# import functions.params
from functions.params import *


#######################################

def get_ic(u0, v0, w0, rho0, T0, u1, v1, w1, rho1, T1):

    for i in range(0, Nzc):
        for j in range(0, Nrf):
            v0[i, j] = 0.0; # v0 = 0.0

    rho0 = np.zeros((Nzc, Nrc));
    T0 = np.ones((Nzc, Nrc));

    for i in range(0, Nzc):
        for j in range(0, Nrc):
            rho0[i, j]     = np.exp(0.5*gamma*M*M*(rc[j]*rc[j]-1.0))

    u1  =   u0;
    v1  =   v0;
    w1  =   w0;

    rho1    = rho0;
    T1      = T0;
    # print("rho0 = \n", rho0);
    return u0, v0, w0, rho0, T0, u1, v1, w1, rho1, T1
#######################################
