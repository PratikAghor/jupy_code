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
        filename = "tau" + n_str
        hf = h5py.File('data/' + filename + '.h5', 'r')
        tau1 = hf.get('tau(t)')
        tau = np.array(tau1)

        # print "\n tau = \n", tau

        ax.plot(tau, z, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = "tau"+ n_str)  # Plot the numerical solution

ax.set_ylabel(r'$z$')  # Set y label
ax.set_xlabel(r'$tau$')  # Set x label
ax.legend()
plt.savefig('tau_t.png')
#######################################
