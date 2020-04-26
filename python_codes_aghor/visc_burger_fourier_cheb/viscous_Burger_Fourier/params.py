import numpy as np  # Import numpy
import math

from numpy import pi, cosh, exp, round, zeros, identity, arange, real, cos, sin, multiply, transpose

#######################################
# define parameters
# Grid and required matrices:
# c stands for computational domain
Lx = 1.0;
xc_min = 0;
xc_max = 2*pi;

# xc = ax + b

a = xc_max - xc_min;
b = xc_min;


nsave = 50; # save after nsave timesteps
Nt = 500; # # of timesteps
nu = 0.01; # viscosity
# theta = 0 => implicit, theta = 1 =>explicit, theta = 0.5 => AM2
theta = 0.5; # weight to the current time-step value for the linear operator.

Nx = 32;

dxc = (xc_max-xc_min)/(Nx-1)

dt = 0.002;

xc = arange(xc_min, xc_max, dxc)
# xc = np.linspace(xc_min, xc_max-dxc, Nx)

print "xc = \n", xc
x = (xc - b*np.ones(Nx-1))/a;

# _padded for plotting
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
