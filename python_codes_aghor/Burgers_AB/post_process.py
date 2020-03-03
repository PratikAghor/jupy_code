# post processing
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import params
from params import *
#######################################
"""
post-process heat_CN data
"""
#######################################
# declare an array that stores the value of u at a point

fig1 = plt.figure(1)  # Create a figure instance
ax1 = fig1.gca()

for n in range(0, AB2_Nt):
    # read data
    if ((n % AB2_nsave) == 0):

        n_str = str(n)
        filenameAB2 = "AB2_u" + n_str
        hfAB2 = h5py.File('data/' + filenameAB2 + '.h5', 'r')
        utAB2 = hfAB2.get('u(t)')
        uAB2 = np.array(utAB2)
        # print u
        ax1.plot(x, uAB2, '-', linewidth=3.0, fillstyle='none', markeredgewidth = 2, label = r'$u$')  # Plot the numerical solution
        ax1.set_xlabel('x')  # Set x label

#ax.legend()
plt.savefig('ut_AB2.png')
#######################################
#######################################
fig2 = plt.figure(2)  # Create a figure instance
ax2 = fig2.gca()
for n in range(0, AB4_Nt):
    # read data
    if ((n % AB4_nsave) == 0):

        n_str = str(n)

        filenameAB4 = "AB4_u" + n_str
        hfAB4 = h5py.File('data/' + filenameAB4 + '.h5', 'r')
        utAB4 = hfAB4.get('u(t)')
        uAB4 = np.array(utAB4)
        ax2.plot(x, uAB4, '-', linewidth=3.0, fillstyle='none', markeredgewidth = 2, label = r'$u$')  # Plot the numerical solution
        ax2.set_xlabel('x')  # Set x label

#ax.legend()
plt.savefig('ut_AB4.png')
#######################################
