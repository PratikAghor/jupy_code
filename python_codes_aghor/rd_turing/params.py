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
nx = 64; ny = 32

xmin = 0.0; xmax = 4*np.pi;
dx = (xmax-xmin)/(float(nx));

ymin = 0.0; ymax = 2*np.pi;
dy = (ymax-ymin)/(float(ny));

x = np.linspace(xmin, xmax-dx, nx); # NOTE: last point xmax not included
y = np.linspace(ymin, ymax-dx, ny); # NOTE: last point xmax not included

x_full = np.linspace(xmin, xmax, nx+1); # for plotting
y_full = np.linspace(ymin, ymax, ny+1); # for plotting

A1x_colloc_pbc, B1x_colloc_pbc = create_matrices_interpolate_1_colloc_pbc(x)
A2x_colloc_pbc, B2x_colloc_pbc = create_matrices_interpolate_2_colloc_pbc(x)

A1y_colloc_pbc, B1y_colloc_pbc = create_matrices_interpolate_1_colloc_pbc(y)
A2y_colloc_pbc, B2y_colloc_pbc = create_matrices_interpolate_2_colloc_pbc(y)

# print('np.shape(B2x_colloc_pbc) =', '\n', np.shape(B2x_colloc_pbc))
# print('np.shape(A2y_colloc_pbc) =', '\n', np.shape(A2y_colloc_pbc))

#######################################
# params
a = 1
b = 5
d = 60

# time stuff
t = 0
dt = 1e-3
Nt = 20000 # # time steps
nsave = 100
#######################################
