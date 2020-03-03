import numpy as np  # Import numpy
import math

import Burger_AB_functions
from Burger_AB_functions import *

#######################################
# define parameters
nu = 0.01; # thermal diffusivity
Nx = 50; # no. of grid pts
u = np.zeros(Nx);
xmin = 0.0; xmax = 1.0; # domain length from [xmin, xmax]
x = np.linspace(xmin, xmax, Nx); # the grid
dx = (xmax - xmin)/(float)(Nx-1.0); # grid size

AB2_dt = 1e-2; # time step for AB2
AB2_Nt = 20; # number of time steps for AB2
AB2_t =  np.linspace(0, AB2_dt*(AB2_Nt-1), AB2_Nt); # time points for AB2

AB2_nsave = AB2_Nt/4; # save data after every nsave steps

AB4_dt = 1e-2; # time step for AB2
AB4_Nt = 20; # number of time steps for AB2
AB4_t =  np.linspace(0, AB4_dt*(AB4_Nt-1), AB4_Nt); # time points for AB2

AB4_nsave = AB4_Nt/4; # save data after every nsave steps

# IC:

# I => global grid, including ghost nodes
uI_0 = np.zeros(Nx+2);

for i in range(0, Nx):
    uI_0[i+1] = np.sin(2.0*(np.pi)*x[i]);
set_periodic_BCs(uI_0);


# declare rhs
rhs0 = np.zeros(Nx+2);
rhs1 = np.zeros(Nx+2);
rhs2 = np.zeros(Nx+2);
rhs3 = np.zeros(Nx+2);
#######################################
