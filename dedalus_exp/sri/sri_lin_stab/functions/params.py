# Computation of Strongly Compressible Rotating Flows, Harada, J. Comp. Phys., 1980;
# Author: Pratik Aghor
#######################################
import numpy as np  # Import numpy
import math


import sys
sys.path.append('/home/aghor/Aghor/UNH/independent/tcf_lin_stab_local')

# for saving figs in this folder
import os
my_path = os.path.abspath(__file__)

from functions import *
from functions.cheb import *
#######################################
#######################################
# parameters
# Re = 100.0; # Reynolds #
E = 1.03e-3; # Ekman #
Pr = 0.97; # Prandtl #
gamma = 5.0/3.0; # gamma = cp/cv
Gamma = 1.0; # aspect ratio - in Hyun and Park = H/L
M = 4.0; # Ma #
oneBygammaMsq = 1.0/(gamma * M * M);
gammaEbyPr = gamma * E * (1.0/Pr);

mu = 1.0; # dimensionles viscosity
epsilon = 3.25e-2; # thermal Rossby number
#######################################
#######################################
"""
define parameters
# Grid and required matrices:
# c stands for computational domain
"""
rmin = 0.3; rmax = 1.0;
zmin = 0.0; zmax = Gamma;

xc_min = -1.0;
xc_max = 1.0;
#######################################
"""
x is the Physical domain [r_i, r_o)
xc is the computational domain [-1, 1]

rc = ar*r + br
zc = az*z + bz
"""

ar = -2.0/(rmax - rmin);
br = 1.0*(rmax + rmin)/(rmax - rmin);
az = 2.0/(zmax - zmin);
bz = -1.0*(zmax + zmin)/(zmax - zmin);


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
Nr = 16; Nz = 16;

D1r, rc = cheb(Nr);
D1z, zc = cheb(Nz);
# print "xc = \n", xc
r = (rc - br*np.ones(Nr+1))/ar;
z = (zc - bz*np.ones(Nz+1))/az;

#######################################
