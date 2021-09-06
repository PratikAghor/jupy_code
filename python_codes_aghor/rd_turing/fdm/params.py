#######################################
# spin-up Hyun and Park JFM, 1992;
# A compact FDM on staggered grid for Navier-Stokes flows, Zhang et. al., Int. J. Numer. Meth. Fluids. 2006;
# Author: Pratik Aghor
import numpy as np
import interpolate
from interpolate import *

#######################################
#######################################
# domain
nx = 64; ny = 16

xmin = 0.0; xmax = 8.0*np.pi;
dx = (xmax-xmin)/(float(nx));

ymin = 0.0; ymax = 2.0*np.pi;
dy = (ymax-ymin)/(float(ny));

Lx = xmax - xmin
Ly = ymax - ymin

x = np.linspace(xmin, xmax-dx, nx); # NOTE: last point xmax not included
y = np.linspace(ymin, ymax-dx, ny); # NOTE: last point xmax not included

x_full = np.linspace(xmin, xmax, nx+1); # for plotting
y_full = np.linspace(ymin, ymax, ny+1); # for plotting

Dx = get_d1x_mat(nx+2, dx)
D2x = get_d2x_mat(nx+2, dx)
Dy = get_d1y_mat(ny+2, dy)
D2y = get_d2y_mat(ny+2, dy)

#######################################
# params
a = 0
b = 1
d = 100

# time stuff
t = 0
dt = 1e-4
Nt = 500000 # # time steps
nsave = 1000
#######################################
