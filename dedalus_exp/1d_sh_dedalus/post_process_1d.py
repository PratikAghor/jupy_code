# post processing
import h5py
import numpy as np
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
from matplotlib.animation import FuncAnimation
from dedalus import public as de

import make_movie
from make_movie import *
#######################################
"""
post-process 1d Swift-Hohenberg equation (sh1d)

dt(u) = ru - [(dx**2 + qc**2)**2] u + v*u**2-g*u**3
"""
#######################################
Nt = 20000;
nsave = Nt/500;
dt = 1.0e-2;

Lx = 200.0
Nx = 512

x_basis = de.Fourier('x', Nx, interval=(-Lx/2, Lx/2), dealias=3/2)
domain = de.Domain([x_basis], np.float64)

x = domain.grid(0)
#######################################
# fig = plt.figure(1)  # Create a figure instance
# ax = fig.gca()
#
# for n in range(0, Nt+1):
#     # read data
#
#     if ((n % nsave) == 0):
#         #print("n =" , n, '\n')
#
#         n_str = str(n)
#         filename = "c" + n_str
#         hf = h5py.File('data/' + filename + '.h5', 'r')
#         ct = hf.get('c(t)')
#         c = np.array(ct)
#
#         ax.plot(x, c, '-', linewidth=2.0, fillstyle='none', markeredgewidth = 2\
#         ,label = "c"+n_str)   # Plot the numerical solution
#
# ax.set_xlabel(r'$x$')  # Set x label
# ax.set_ylabel(r'$c$')  # Set y label
# ax.legend()
# plt.savefig('ct_fourier.png')
#######################################
make_dedalus_movie(Nt, nsave, Lx, Nx, [-1, 1], varname="u", name="shper_2", show=0)
#######################################
