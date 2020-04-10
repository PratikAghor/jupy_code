import cheb
from cheb import *

from numpy import pi, cos, sin

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
###########################################
'''
Test cheb differentiation
author: Pratik Aghor
'''
###########################################
N = 50;

D, xc = cheb(N);

# print "D = \n ", D
u = sin(pi*xc); # Initial condition

u_1x_exact = pi* cos(pi*xc)

u_1x_cheb = np.matmul(D, u)


###########################################
fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()

ax.plot(xc, u_1x_cheb, 'o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "cheb")  # Plot the numerical solution
ax.plot(xc, u_1x_exact, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "exact")  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$u_{1x}$')  # Set y label
ax.legend()
plt.savefig('test_cheb.png')
###########################################
