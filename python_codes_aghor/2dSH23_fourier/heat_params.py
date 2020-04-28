import numpy as np  # Import numpy
import math

from numpy import pi, cosh, exp, round, zeros, identity, arange, real, cos, sin,\
 multiply, transpose, mean, sqrt

#######################################
# define parameters
# Grid and required matrices:
# c stands for computational domain
Lx = 5; Ly = 5;
xc_min = 0.;
xc_max = 2.0*pi;
yc_min = 0.;
yc_max = 2.0*pi;

# xc = a1x + b1 ; yc = a2x + b2

a1 = (xc_max - xc_min)*1.0/Lx;
b1 = xc_min;

a2 = (yc_max - yc_min)*1.0/Ly;
b2 = yc_min;


nsave = 100; # save after nsave timesteps
Nt = 1000; # 40000; # # of timesteps

nu = 2.2; # 1.6; # N(u) = nu u^2 - u^3
mu = 0.5; # 0.31;
kappa = 0.1; # for the heat equation test
# theta = 0 => implicit, theta = 1 =>explicit, theta = 0.5 => AM2
theta = 0.5; # weight to the current time-step value for the linear operator.

Nx = 32; Ny = 32;

Nxc = Nx-1; Nyc = Ny - 1;

dxc = (xc_max-xc_min)/(Nxc)
dyc = (yc_max-yc_min)/(Nyc)

dt = 0.1;

xc = arange(xc_min, xc_max, dxc)
yc = arange(yc_min, yc_max, dyc)
# print ("xc = \n", xc)
x = (xc - b1*np.ones(Nxc))/a1;
y = (yc - b2*np.ones(Nxc))/a2;

# _padded for plotting
xc_padded = arange(xc_min, xc_max+dxc, dxc)
yc_padded = arange(yc_min, yc_max+dyc, dyc)

x_padded = (xc_padded - b1*np.ones(Nx))/a1;
y_padded = (yc_padded - b2*np.ones(Ny))/a2;


kx = zeros(Nxc);
kx[0:Nx/2] = arange(0, Nx/2); kx[Nx/2:] = arange(-Nx/2+1, 0, 1);

ky = zeros(Nyc);
ky[0:Ny/2] = arange(0, Ny/2); ky[Ny/2:] = arange(-Ny/2+1, 0, 1);

# print "kx= \n", kx
# print "ky= \n", ky
#######################################
# Initial condition u0
# u0 = zeros((Nyc, Nxc)); hexa = zeros((Nyc, Nxc))
# for i in range(0, Nyc):
#     for j in range(0, Nxc):
#         u0[i, j] = (1.0 + sin(xc[j])*sin(yc[i]));
# # print "u(x, y) = sin(xc[j])*sin(yc[i]) \n"

# u = 4*sech(sqrt(xx.^2+yy.^2)/10).*u;
# for i in range(0, Nyc):
#     for j in range(0, Nxc):
#         u0[i, j] = (2.0/ cosh(sqrt((xc[j]+20)**2 + (yc[i]+20)**2)/10.0) )
#         hexa[i, j] = (cos(xc[j])+ cos((xc[j]+sqrt(3)*yc[i])/2)+ cos((xc[j]-sqrt(3)*yc[i])/2))/3
#         u0[i, j] = u0[i, j] + (2.0/cosh(sqrt(xc[j]**2+yc[i]**2)/5))*hexa[i, j];
#######################################
