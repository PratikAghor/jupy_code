import numpy as np  # Import numpy
import math

import cheb
from cheb import *
#######################################
#######################################
"""
define parameters
# Grid and required matrices:
# c stands for computational domain
"""
Lx = 1.0;
xc_min = -1.0;
xc_max = 1.0;
#######################################
"""
x is the Physical domain [0, 1)
xc is the computational domain [-1, 1]

xc = ax + b
"""

a = xc_max - xc_min;
b = xc_min;


#######################################
"""
nsave: save data after nsave timesteps
Nt   : # of timesteps
nu   : viscosity
dt   : time-step
"""

Nt = 100;
nsave = Nt/5;
nu = 0.01;
dt = 0.01;
#######################################
"""
Nx = # of grid points (kind of)

For given Nx, the Gauss-Lobatto grid will be consist of (Nx+1) points.

D is the Chebyeshev differentiation matrix defined on xc in [-1, 1].

Finally, we can get back to the physical grid by inverting the relation between
x and xc to have: x = (xc - b)/ a .
"""
Nx = 127;

D, xc = cheb(Nx);
# print "xc = \n", xc
x = (xc - b*np.ones(Nx+1))/a;

# print "x = \n", x
#######################################
