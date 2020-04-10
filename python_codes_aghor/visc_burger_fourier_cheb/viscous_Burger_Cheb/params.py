import numpy as np  # Import numpy
import math

import cheb
from cheb import *
#######################################
# define parameters
# Grid and required matrices:
# c stands for computational domain
Lx = 1.0;
xc_min = -1.0;
xc_max = 1.0;

# xc = ax + b

a = xc_max - xc_min;
b = xc_min;


nsave = 50; # save after nsave timesteps
Nt = 500; # # of timesteps
nu = 0.01; # viscosity

Nx = 63;

dt = 0.002;

D, xc = cheb(Nx);
# print "xc = \n", xc
x = (xc - b*np.ones(Nx+1))/a;

# print "x = \n", x
#######################################
