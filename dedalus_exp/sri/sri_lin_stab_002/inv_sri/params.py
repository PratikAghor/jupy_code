# Computation of Strongly Compressible Rotating Flows, Harada, J. Comp. Phys., 1980;
# Author: Pratik Aghor
#######################################
import numpy as np  # Import numpy
import math

# for saving figs in this folder
import os
my_path = os.path.abspath(__file__)

from cheb import *
#######################################
# geometric parameters
Gamma = 10.0; # aspect ratio - in Hyun and Park = H/L
eta = 0.9; # eta = r_i/r_o, radius ratio
mu = 0.4; # Omega_o/Omega_i
Re = 100;
Fr = 1.0;
oneByFrsq = 1.0/(Fr**2)
prefactor = (eta/(1.0-eta))*(1.0/Re)
#######################################
"""
define parameters
# Grid and required matrices:
# c stands for computational domain
"""
rmin = 1.0; rmax = rmin/eta;
xc_min = -1.0;
xc_max = 1.0;
#######################################
"""
x is the Physical domain [r_i, r_o)
xc is the computational domain [-1, 1]

rc = ar*r + br
zc = az*z + bz
"""

ar = -2.0/(rmax - rmin); # - sign makes sure that rc = -1 corresponds to the inner boundary
br = 1.0*(rmax + rmin)/(rmax - rmin);
#######################################

#######################################
"""
Nr = # of grid points (kind of)
Nk = # of values of k for which the calculation will be done

For given Nx, the Gauss-Lobatto grid will be consist of (Nx+1) points.

D is the Chebyeshev differentiation matrix defined on xc in [-1, 1].

Finally, we can get back to the physical grid by inverting the relation between
x and xc to have: x = (xc - b)/ a .
"""
Nr = 16
Nk = 32

D1rc, rc = cheb(Nr);

r = (rc - br*np.ones(Nr+1))/ar
print("r = \n", r)
oneByr = 1.0/r

Id = np.identity(Nr+1)
D1r = ar*D1rc
D2r = np.matmul(D1r, D1r)
#######################################

#######################################
# other parameters
m = 1

kvec = np.zeros(Nk);
kvec[0:Nk] = (2.0*pi/Gamma)*arange(0, Nk)
print("kvec = \n", kvec)
#######################################
# base state description
A = (mu - eta**2)/(1.0 - eta**2)
B = (eta**2)*(1.0-mu)/((1.0 + eta)*(1.0-eta)**3)
Z = 2.0*A
print("A, B, Z = ", A, B, Z)
Omega = A + B/(r*r)

print("Omega = \n", Omega)
#######################################
