import numpy as np  # Import numpy
import math

import heat_CN_functions
from heat_CN_functions import *

#######################################
# define parameters
nu = 0.01; # thermal diffusivity
Nx = 51; # no. of grid pts
xmin = 0.0; xmax = 1.0; # domain length from [xmin, xmax]
x = np.linspace(xmin, xmax, Nx); # the grid
dx = (xmax - xmin)/(float)(Nx-1.0); # grid size

dt = 0.1; # time step
Nt = 100; # number of time steps
t =  np.linspace(0, dt*(Nt-1), Nt); # time points
nsave = 5; # save data after every nsave steps


alpha = (nu*dt)/(float)(2.0*dx*dx); # alpha - convenience variable

u_l = 0.0; # value of u at the left boundary
u_r = 0.0; # value of u at the right boundary

# IC:
u_0 = np.sin(4.0*(np.pi)*x);
set_Dirichlet_BCs(u_0, u_l, u_r);
#######################################
