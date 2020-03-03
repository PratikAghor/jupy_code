import numpy as np  # Import numpy
import math
import h5py # for saving data
from numpy.linalg import *

import Burger_AB_functions
from Burger_AB_functions import *

import params
from params import *
#######################################
'''
Code to integrate 1d Burger's equation (u_t + uu_x = 0) with periodic boundaries
using finite difference method. AB2 and AB4 were used for time marching.
Author: Pratik Aghor
'''
#######################################
# IC:
# u_0 defined in params
# print "np.shape(uI_0) = " '\t', np.shape(uI_0)
# print "x = ", '\n', x
# print "uI_0 = ", '\n', uI_0
uI_1 = uI_0;
#######################################
# time marching

for n in range(0, AB2_Nt+1):
    # impelement BC's
    set_periodic_BCs(uI_0);
    set_periodic_BCs(uI_1);

    # get rhs
    rhs0 = get_Burger_rhs(uI_0, dx);
    rhs1 = get_Burger_rhs(uI_1, dx);
    # print "rhs0 = ", '\n', rhs0

    # AB2 time marching
    uI_2 = uI_1 + 0.5*AB2_dt*(3.0*rhs1 - rhs0);

    # save data after every nsave steps
    if ((n % AB4_nsave) == 0):
        for i in range(0, Nx):
            u[i] = uI_1[i+1]
        n_str = str(n)
        filename = "AB2_u" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('u(t)', data=u)
        hf.close()

    # update u_nplus1 to u_n
    uI_0 = uI_1;
    uI_1 = uI_2;
#######################################
