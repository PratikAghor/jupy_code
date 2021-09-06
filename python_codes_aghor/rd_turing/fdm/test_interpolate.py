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
n = 64; # # of center points

xmin = 0.0; xmax = 2.0*np.pi;
dx = (xmax-xmin)/(float(n));

x = np.linspace(xmin, xmax-dx, n); # NOTE: last point xmax not included
x_full = np.linspace(xmin, xmax, n+1); # for plotting
print('x =', '\n', x)
print("dx =",  dx)
#######################################
# test interpolate_0

f = np.sin(x)
f_padded = np.zeros(n+2) # ghost nodes included
f_padded[1:n+1] = f
# periodic bc (pbc)
f_padded[0] = f_padded[n]
f_padded[n+1] = f_padded[1]
print("f = ", f)
print("f_padded = ", f_padded)

f_1x_exact = np.cos(x)
f_2x_exact = -np.sin(x)

D = get_d1x_mat(n+2, dx)
D2 = get_d2x_mat(n+2, dx)

# print("D = \n", D)
# print("D2 = \n", D2)

f_1x = D @ f_padded
f_2x = D2 @ f_padded

print("f_1x = \n", f_1x)
#######################################
#######################################
fig1 = plt.figure(1)  # Create a figure instance
ax1 = fig1.gca()  # projection='3d' to Get current axes in 3D projection
ax1.plot(x, f, '-o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f$' )  # Plot the solution
ax1.plot(x, f_1x[1:n+1], '^', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f_{1x}$ numerical' )
ax1.plot(x, f_1x_exact, '--', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f_{1x}$ exact' )  # Plot the solution
ax1.set_xlabel('x')  # Set x label
ax1.legend()
plt.savefig('test_interpolate_1.png')
#######################################
#######################################
fig2 = plt.figure(2)  # Create a figure instance
ax2 = fig2.gca()  # projection='3d' to Get current axes in 3D projection
ax2.plot(x, f, '-o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f$' )  # Plot the solution
ax2.plot(x, f_2x[1:n+1], '^', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f_{2x}$ numerical' )
ax2.plot(x, f_2x_exact, '--', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$f_{2x}$ exact' )  # Plot the solution
ax2.set_xlabel('x')  # Set x label
ax2.legend()
plt.savefig('test_interpolate_2.png')
#######################################
