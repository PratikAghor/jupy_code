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
post-process AM2AB2 2dheat/2dSH data
"""
#######################################

fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()

for n in range(0, Nt+1):

    t = 0 + dt*n;
    t_str = str(t)
    # read data
    if ((n % nsave) == 0):

        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('heat_data/' + filename + '.h5', 'r')
        ut = hf.get('u(t)')
        u_padded = np.array(ut)

        ax.plot(x_padded, u_padded[(Ny/2)-1, :], '-', linewidth=2.0, fillstyle='none', \
        markeredgewidth = 2, label = "u_" + t_str)  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$u$')  # Set y label
ax.legend()
plt.savefig('ut_heat_fourier.png')
#######################################
