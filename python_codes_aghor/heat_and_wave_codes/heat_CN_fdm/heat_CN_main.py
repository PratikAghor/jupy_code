import numpy as np  # Import numpy
import math
import h5py # for saving data
from numpy.linalg import *

import heat_CN_functions
from heat_CN_functions import *

import params
from params import *
#######################################
'''
Code to integrate 1d heat equation (u_t = nu*u_xx) with Dirichlet boundaries
using finite difference method. Crank-Nicholson CN was used for time marching.
Author: Pratik Aghor
'''
#######################################
# print "alpha = ", alpha
# print "1+ 2*alpha = ", 1.0+2.0*alpha
# print "1- 2*alpha = ", 1.0-2.0*alpha

# IC:
# u_0 defined in params
u_n = u_0;

# get LHS and RHS matrices
A = build_A(Nx, alpha);
B = build_B(Nx, alpha);

# print "A = ", '\n', A
# print "B = ", '\n', B

rhs0 = np.matmul(B, u_n);
# print "rhs0", '\n', np.transpose(rhs0)

#######################################
# time marching

for n in range(0, Nt+1):
    # impelement BC's
    set_Dirichlet_BCs(u_n, u_l, u_r);

    rhs = np.matmul(B, u_n); # build rhs vector
    u_nplus1 = qrsolve(A, rhs);

    # save data after every nsave steps
    if ((n % nsave) == 0):
        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('u(t)', data=u_n)
        hf.close()

    # update u_nplus1 to u_n
    u_n = u_nplus1;
#######################################
