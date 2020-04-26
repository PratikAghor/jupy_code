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

# Initial condition
u = sin(xc); 

u_padded = zeros(Nx);

alpha = 1j*k;
# print "alpha=", alpha

# Writing some operators for convenience
# # D = d/d(xc) operator in Fourier space
# L = linear operator nu*D^2 in Fourier space, diagonal matrix!
# G = -1/2 D operator in Fourier space

D = a*(alpha); 
L = nu*D*D;                                     
G = -0.5*D;                                     

# A     = (I - dt/2*L)
A_inv = (np.ones(Nx-1) - (1.0-theta)*dt*L)**(-1);

B     = np.ones(Nx-1) + theta*dt*L;

# writing nonlinear term 
#-u u_x (spectral), notation N1 = N^n = N(u(n dt))
#notation N0 = N^{n-1} = N(u((n-1) dt))
N1  = G*fft(u*u);                               
N0 = N1;                                        

# transform u to spectral
v  = fft(u);                                   
#######################################
# time marching!
for n in range (0, Nt):
    N0 = N1;
    N1 = G*fft(u*u); # compute N0 = -u u_x


    v = A_inv*((B*v) + 1.5*dt*N1 - 0.5*dt*N0);

    u = real(ifft( v ))
	
	# write u_padded with the last term = u[0] on Nx grid. for plotting purposes 
    u_padded[0:Nx-1] = u[0:Nx-1]; 
    u_padded[Nx-1] = u[0]; 
    # save data after every nsave steps
    if ((n % nsave) == 0):
        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('u(t)', data=u_padded)
        hf.close()
