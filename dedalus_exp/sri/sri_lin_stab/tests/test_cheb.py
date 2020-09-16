from numpy import pi, cos, sin

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

import sys
sys.path.append('/home/aghor/Aghor/UNH/independent/jupy_codes_local/dedalus_exp/sri/sri_lin_stab')

# for saving figs in this folder
import os
my_path = os.path.abspath(__file__)

from functions import *
from functions.cheb import *

###########################################
'''
Test cheb differentiation
author: Pratik Aghor
'''
###########################################
N = 50;

D, xc = cheb(N);

# print "D = \n ", D
n = 3
u = sin(n * pi * xc); # Initial condition
u_1x_exact = n * pi * cos(n * pi * xc)
u_1x_cheb = np.matmul(D, u)

f = np.exp(xc)
f_1x_exact = f
f_1x_cheb = np.matmul(D, f)
###########################################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()

ax.plot(xc, u_1x_cheb, 'o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "cheb")  # Plot the numerical solution
ax.plot(xc, u_1x_exact, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "exact")  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$u_{1x}$')  # Set y label
ax.legend()
plt.savefig(my_path + 'test_cheb_1.png')
###########################################
fig = plt.figure(2)  # Create a figure instance
ax = fig.gca()

ax.plot(xc, f_1x_cheb, 'o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "cheb")  # Plot the numerical solution
ax.plot(xc, f_1x_exact, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "exact")  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$f_{1x}$')  # Set y label
ax.legend()
plt.savefig(my_path + 'test_cheb_2.png')
###########################################
