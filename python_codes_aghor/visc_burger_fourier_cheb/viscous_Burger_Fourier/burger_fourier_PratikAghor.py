from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from numpy import pi, cosh, exp, round, zeros, identity, arange, real, cos, sin, multiply, transpose
from numpy.fft import fft,ifft
from matplotlib.pyplot import figure

import time
import numpy as np  # Import numpy
import math
from numpy.linalg import inv
import h5py # for saving data

import params
from params import *

###########################################
'''
Code to integrate 1d Burger's equation (u_t + u*u_x = nu*u_xx) with periodic boundaries using Fourier spectral method.
For time integration, the second-order (two step) Adams-Moulton method for the linear term
and the second-order (two-step) Adams-Bashforth method for the nonlinear term are used (Henceforth called AM2AB2).

The code is written on the transformed co-ordinates xc = [0, 2\pi]
Physical coordinates [0, 1]

- Pratik Aghor
'''
###########################################
# print "x=", x

# Initial condition
u0 = sin(xc);

u0_padded = zeros(Nx); u1_padded = zeros(Nx);

alpha = 1j*k;
###########################################
"""
Writing some operators for convenience
D = d/d(xc) operator in Fourier space
L = linear operator nu*D^2 in Fourier space, diagonal matrix!
G = -1/2 D operator in Fourier space
A     = (I - dt/2*L)
B 	= (I + dt/2*L)

NOTE: Since L is a diagonal matrix, we store A and B as 1d vectors.

They are actually stand-ins for an (Nx-1) x (Nx-1) diagonal matrices.

Then the matrix mult B*x becomes element-wise mult for the stored 1d vector.
"""

D = a*(alpha);
L = nu*D*D;
G = -0.5*D;

A_inv = (np.ones(Nx-1) - (1.0-theta)*dt*L)**(-1);

B     = np.ones(Nx-1) + theta*dt*L;
###########################################
"""
writing nonlinear term
-u u_x (spectral), notation N1 = N^n = N(u(n dt))
notation N0 = N^{n-1} = N(u((n-1) dt))

We initialize with u1 = u0, then the first step in AM2AB2 algorithm
changes into a CN-Euler step, and from second step, it is AM2AB2 time-marching
"""
u1  = u0;
N0  = G*fft(u0*u0);
N1 = N0;
#######################################
"""
time marching!

Given u0 and u1, find u2 and update

v0 and v1 are u0 an u1 in spectral space v0  = fft(u0); v1  = fft(u1);

"""

t1 = time.time()
for n in range (0, Nt+1):
    v0  = fft(u0); v1  = fft(u1);
    N0  = G*fft(u0*u0); N1 = G*fft(u1*u1);

    v2 = A_inv*((B*v1) + 1.5*dt*N1 - 0.5*dt*N0);
    u2 = real(ifft( v2 ));

    # write u_padded with the last term = u[0] on Nx grid. for plotting purposes
    u1_padded[0:Nx-1] = u1[0:Nx-1];
    u1_padded[Nx-1] = u1[0];

    # save data after every nsave steps
    if ((n % nsave) == 0):
        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('u(t)', data=u1_padded)
        hf.close()

    # update variables for next time-step
    u0 = u1;
    u1 = u2;

t2 = time.time()
print "time required for Fourier Spectral = ", t2 - t1, '\n'
