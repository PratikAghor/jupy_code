from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from numpy import pi, cosh, exp, round, zeros, identity, arange, real, cos, sin, multiply, transpose
from numpy.fft import fft,ifft
from matplotlib.pyplot import figure

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
u = sin(xc); # Initial condition

alpha = 1j*kx; # real wavenumbers:    exp(alpha*x)
# print "alpha=", alpha

D = a*(alpha); # D = d/d(xc) operator in Fourier space
# print "D =", '\n', D;

L = nu*D*D;                                     # linear operator nu*D^2 in Fourier space, diagonal matrix!

G = -0.5*D;                                     # -1/2 D operator in Fourier space

# A     = (np.ones(Nx) - dt2*L)
# print "A = ", '\n', A
A_inv = (np.ones(Nx) - (1.0-theta)*dt*L)**(-1);
# print "A_inv = ", '\n', A_inv
B     = np.ones(Nx) + theta*dt*L;

N1  = G*fft(u*u);                               #-u u_x (spectral), notation N1 = N^n = N(u(n dt))
N0 = N1;                                        #notation N0 = N^{n-1} = N(u((n-1) dt))
v  = fft(u);                                    # transform u to spectral
#######################################
# time marching!
for n in range (0, Nt):
    N0 = N1;
    N1 = G*fft(u*u); # compute N0 = -u u_x


    v = A_inv*((B*v) + 1.5*dt*N1 - 0.5*dt*N0);

    u = real(ifft( v ))
    # save data after every nsave steps
    if ((n % nsave) == 0):
        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('u(t)', data=u)
        hf.close()
