import numpy as np  # Import numpy
import math

from numpy import pi, cosh, exp, round, zeros, identity, arange, real, cos, sin, multiply, transpose

#######################################
"""
define parameters
# Grid and required matrices:
# c stands for computational domain
"""
Lx = 1.0;
xc_min = 0;
xc_max = 2*pi;
#######################################
"""
x is the Physical domain [0, 1)
xc is the computational domain [0, 2\pi)

NOTE: We do not include the last point here, because periodic BC's are assumed.
Also, fft assumes periodic BC's and it would be erroneous to do the calculation
inclding the last point.

While saving and plotting, I pad up the solution as well as get the x_padded, in
order to get the solution on the whole domain [0, 1]

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

theta : weight to the current time-step value for the linear operator.
if theta = 0 => implicit, if theta = 1 => explicit, theta = 0.5 => AM2
here, theta = 0.5 is used.
"""

Nt = 100;
nsave = Nt/5;
nu = 0.01;
theta = 0.5;
dt = 0.01;
#######################################
"""
Nx : # of grid points

dxc = dx on the computational grid [0, 2\pi)

_padded is used for saving and plotting

kx is the one with uneven distribution of wave-numbers,
for example, for Nx = 8,

kx = [0, 1 , 2, 3, 4, -3, -2, -1]
To make it balanced (see blog for more discussion on this),
we zero-out the highest mode:

Hence kx = [0, 1, 2, 3, 0, -3, -2, -1]

However, since we have periodic BC's assumed and we are neglecting the last grid
point, we want an (Nx-1) sized kvector

Therefore, k = [0, 1, 2, 3, -3, -2, -1].
"""
Nx = 128;

dxc = (xc_max-xc_min)/(Nx-1)

xc = arange(xc_min, xc_max, dxc)

x = (xc - b*np.ones(Nx-1))/a;


xc_padded = arange(xc_min, xc_max+dxc, dxc)
x_padded = (xc_padded - b*np.ones(Nx))/a;

kx = zeros(Nx);
kx[0:Nx/2] = arange(0,Nx/2); kx[Nx/2+1:] = arange(-Nx/2+1,0,1);
# print "kx=", kx

k = zeros(Nx-1);
for i in range(0, Nx-1):
	if i < Nx/2:
		k[i] = kx[i];
	else:
		k[i] = kx[i+1];

#######################################
