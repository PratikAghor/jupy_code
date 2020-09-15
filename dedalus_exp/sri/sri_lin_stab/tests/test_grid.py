from numpy import pi, cos, arange, ones

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

import sys
sys.path.append('/home/aghor/Aghor/UNH/independent/tcf_lin_stab_local')

# for saving figs in this folder
import os
my_path = os.path.abspath(__file__)

from functions import *

###########################################
'''
Test cheb grid
author: Pratik Aghor
'''
###########################################
# print("we make r go from ra to rb and z from Gamma to 0")
# print("so that increasing i goes from Gamma to 0")
# print("increasing j goes from ra to rb")
# print("to keep the matrices aligned with the physical picture.\n")
# print("\n r = (should be between ra and rb)\n", r)
# # print("\n z = (should be between Gamma and 0)\n", z)
#
# print("\n rc = (should be between 1 and -1)\n", rc)
# print("\n zc = (should be between 1 and -1)\n", zc)
###########################################
N = 5

xc  = cos(pi*arange(0,N+1)/N) # xc = Chebyshev-Gauss-Lobatto grid
xg  = cos(pi*arange(1, (2*N), 2)/(2.0*(N)))  #   xg = Gauss grid

print("xc = ", xc)
print("xg = ", xg)
###########################################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()

ax.plot(xc, np.zeros(N+1), 'o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "Cheb-Gauss-Lobatto")  # Plot the numerical solution
ax.plot(xg, np.zeros(N), '^', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "Cheb-Gauss")  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
# ax.set_ylabel(r'$u_{1x}$')  # Set y label
ax.legend()
plt.savefig('tests/test_gauss_grid.png')
###########################################
