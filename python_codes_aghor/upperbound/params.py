import numpy as np  # Import numpy
from numpy import pi, arange, cos, sin, zeros
import math

import cheb
from cheb import *
#######################################
# define parameters
#######################################
# Grid and required matrices:
# c stands for computational domain
Lz = 1.0; Lx = 2.0;
zc_min = -1.0;
zc_max = 1.0;

# zc = az + b
a = zc_max - zc_min;
b = zc_min;

nsave = 5000; # save after nsave timesteps
Nt = 10000; # # of timesteps
dt = 1e-4;

Nz = 64;
Nx = 12;

D1, zc = cheb(Nz);
D2 = np.matmul(D1, D1)
# print "xc = \n", xc
z = (zc - b*np.ones(Nz+1))/a;

# print "z = \n", z
#######################################
Ra = 50.0; # Rayleigh #
ra = 2.0 * Ra;


kvec = zeros(Nx);
kvec[0:Nx] = (2.0*pi/Lx)*arange(1,Nx+1)
#######################################
