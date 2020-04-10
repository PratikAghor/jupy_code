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
post-process AM2AB2 Burger's data
"""
#######################################

fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()

for n in range(0, Nt):
    # read data
    if ((n % nsave) == 0):

        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('data/' + filename + '.h5', 'r')
        ut = hf.get('u(t)')
        u = np.array(ut)

        ax.plot(x, u, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2)  # Plot the numerical solution

ax.set_xlabel(r'$x$')  # Set x label
ax.set_ylabel(r'$u$')  # Set y label
# ax.legend()
plt.savefig('ut_fourier.png')
#######################################
