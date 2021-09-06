#######################################
# spin-up Hyun and Park JFM, 1992;
# A compact FDM on staggered grid for Navier-Stokes flows, Zhang et. al., Int. J. Numer. Meth. Fluids. 2006;
# Author: Pratik Aghor
import numpy as np
import interpolate
from interpolate import *

#Import plotting functions:
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
#######################################
n = 32; # # of center points

xmin = 0.0; xmax = 2.0*np.pi;
dx = (xmax-xmin)/(float(n));

x = np.linspace(xmin, xmax-dx, n); # NOTE: last point xmax not included
print('x =', '\n', x)

A1_colloc_pbc, B1_colloc_pbc = create_matrices_interpolate_1_colloc_pbc(x)
A2_colloc_pbc, B2_colloc_pbc = create_matrices_interpolate_2_colloc_pbc(x)

# print("shape(A1_colloc_pbc) = ", np.shape(A1_colloc_pbc))
# print('A2_colloc_pbc =', '\n', A2_colloc_pbc)
# print('B2_colloc_pbc =', '\n', B2_colloc_pbc)
#######################################
# test interpolate_0

f = np.sin(x)
f_1x_exact = np.cos(x)
f_2x_exact = -np.sin(x)

f_1_colloc = interpolate_colloc(f, x, A1_colloc_pbc, B1_colloc_pbc)
f_2_colloc = interpolate_colloc(f, x, A2_colloc_pbc, B2_colloc_pbc)

#######################################
#######################################
fig1 = plt.figure(1)  # Create a figure instance
ax1 = fig1.gca()  # projection='3d' to Get current axes in 3D projection
ax1.plot(x, f, '-o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f$' )  # Plot the solution
ax1.plot(x, f_1_colloc, '^', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f_{1x}$ numerical' )
ax1.plot(x, f_1x_exact, '--', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f_{1x}$ exact' )  # Plot the solution
ax1.set_xlabel('x')  # Set x label
ax1.legend()
plt.savefig('test_interpolate_1_colloc_pbc.png')
#######################################
#######################################
fig2 = plt.figure(2)  # Create a figure instance
ax2 = fig2.gca()  # projection='3d' to Get current axes in 3D projection
ax2.plot(x, f, '-o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f$' )  # Plot the solution
ax2.plot(x, f_2_colloc, '^', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f_{2x}$ numerical' )
ax2.plot(x, f_2x_exact, '--', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f_{2x}$ exact' )  # Plot the solution
ax2.set_xlabel('x')  # Set x label
ax2.legend()
plt.savefig('test_interpolate_2_colloc.png')
#######################################
