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

# fig = plt.figure(1)  # Create a figure instance
# ax = fig.gca()
#
# for n in range(0, Nt+1):
#     # read data
#     if ((n % nsave) == 0):
#
#         n_str = str(n)
#         filename = "u" + n_str
#         hf = h5py.File('data/' + filename + '.h5', 'r')
#         ut = hf.get('u(t)')
#         u_padded = np.array(ut)
#
#         ax.plot(x_padded, u_padded[(Ny/2)-1, :], '-', linewidth=2.0, fillstyle='none', \
#         markeredgewidth = 2)  # Plot the numerical solution
#
# ax.set_xlabel(r'$x$')  # Set x label
# ax.set_ylabel(r'$u$')  # Set y label
# # ax.legend()
# plt.savefig('ut_fourier.png')
#######################################
#######################################
fig = plt.figure(1)  # Create a figure instance
ax = fig.add_subplot(111)

for n in range(0, Nt+1):
    # read data
    if ((n % nsave) == 0):

        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('data/' + filename + '.h5', 'r')
        ut = hf.get('u(t)')
        u_padded = np.array(ut)

        plt.imshow(u_padded)  # Plot the numerical solution
        ax.set_aspect('equal')

        # cax = fig.add_axes([0, Nx, 0, Ny])
        # cax.get_xaxis().set_visible(False)
        # cax.get_yaxis().set_visible(False)
        # cax.patch.set_alpha(0)
        # cax.set_frame_on(False)        # ax.legend()
        plt.savefig('snapshots/u' + n_str + '.png')
#######################################
