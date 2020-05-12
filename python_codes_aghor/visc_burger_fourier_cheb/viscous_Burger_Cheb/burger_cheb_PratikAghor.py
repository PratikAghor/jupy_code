from matplotlib.collections import LineCollection
from numpy import pi, arange, cos, sin
from numpy.fft import fft,ifft
from matplotlib.pyplot import figure

import numpy as np  # Import numpy
import time
import math
from numpy.linalg import inv
import h5py # for saving data

import params
from params import *

import cheb
from cheb import *

import qrsolve
from qrsolve import *
###########################################
'''
Code to integrate 1d Burger's equation (u_t + u*u_x = nu*u_xx) with Dirichlet boundaries using Chebyeshev spectral method.
For time integration, the second-order (two step) Adams-Moulton method for the linear term
and the second-order (two-step) Adams-Bashforth method for the nonlinear term are used (Henceforth called AM2AB2).

The code is written on the transformed co-ordinates xc = [-1, 1]
Physical coordinates [0, 1]

- Pratik Aghor
'''
###########################################
# print "x=", x
u0 = -sin(pi*xc); # Initial condition
u1 = u0;


D2 =  np.matmul(D, D) # d^2()/d(xc)^2 on the cheb grid
L = (a**2)*nu*D2;
A = np.identity(Nx+1) - 0.5*dt*L;

# BC
A[0, 0] = 1.0; A[0, 1:-1] = 0.0;
A[Nx, 0:Nx-1] = 0.0; A[Nx, Nx] = 1.0;
#######################################
# time marching!

t1 = time.time()
for n in range (0, Nt):

    u0_sq = np.multiply(u0, u0);
    N0   = -0.5*a*np.matmul(D, u0_sq); # N0 = -(a/2)*(d(u0^2)/dxc)

    u1_sq = np.multiply(u1, u1);
    N1   = -0.5*a*np.matmul(D, u1_sq); # N1 = -(a/2)*(d(u1^2)/dxc)

    rhs = np.matmul( (np.identity(Nx+1) + 0.5*dt*L), u1);
    rhs = rhs + 1.5*dt*N1 - 0.5*dt*N0;

    # bc
    rhs[0] = 0.0; rhs[Nx] = 0.0;

    # solve A u2 = rhs
    u2 = solve(A, rhs);

    # update
    u0 = u1;
    u1 = u2;

    # save data after every nsave steps
    if ((n % nsave) == 0):
        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('u(t)', data=u1)
        hf.close()
#######################################
t2 = time.time()

print "time for Cheb spectral = ", t2 - t1, '\n'

#######################################
