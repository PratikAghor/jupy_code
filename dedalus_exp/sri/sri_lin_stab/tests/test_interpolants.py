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
from functions.interpolants import *

###########################################
'''
Test interpolants
g: Chebyshev-Gauss grid cos((2j+1)pi/2(N+1))
c: Chebyshev-Gauss-Lobatto grid cos(j*pi/N)
j = {0, 1, 2, ..., N}
author: Pratik Aghor
'''
###########################################
N = 63;

G2GL = g2gl(N)
GL2G = gl2g(N)
E = D_g2gl(N)
D, xc = cheb(N);
xg  = cos(pi*arange(1, (2*N), 2)/(2.0*(N)))  # xg = Chebyshev-Gauss grid

# print("G2GL = \n", G2GL)
# print("GL2G = \n", GL2G)

# print "D = \n ", D
n = 3
u = sin(n * pi * xc); # Initial condition
u_1x_exact = n * pi * cos(n * pi * xc)
ug = np.matmul(GL2G, u)

u_g2gl_1x = np.matmul(E, ug)

# print("ug = \n", ug)
# print("len(ug)= \n", len(ug))
# print("ug[0:-1] = \n", ug[0:-1])
###########################################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()

ax.plot(xg, ug[0:-1], 'o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "GL2G")  # Plot the numerical solution
ax.plot(xc, u, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "exact")  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$u$')  # Set y label
ax.legend()
plt.savefig('tests/test_gl2g_interpolant.png')
###########################################
uc = np.matmul(G2GL, ug)
fig = plt.figure(2)  # Create a figure instance
ax = fig.gca()

ax.plot(xc, uc, 'o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "G2GL")  # Plot the numerical solution
ax.plot(xc, u, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "exact")  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$u$')  # Set y label
ax.legend()
plt.savefig('tests/test_g2gl_interpolant.png')
###########################################
###########################################
fig = plt.figure(3)  # Create a figure instance
ax = fig.gca()

ax.plot(xc, u_g2gl_1x, 'o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "D_G2GL")  # Plot the numerical solution
ax.plot(xc, u_1x_exact, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "exact")  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$u$')  # Set y label
ax.legend()
plt.savefig('tests/test_D_g2gl_interpolant.png')
###########################################

###########################################
print("done! Check figures in the tests directory.")
###########################################
