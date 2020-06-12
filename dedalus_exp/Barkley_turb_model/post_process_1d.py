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
post-process Barkley's puffs and slugs model

q = turbulence intensity
u = max velocity in the pipe direction (high for laminar, small for transition)
(as for puffs and slugs, the profile gets blunted)

dt(q) + (u-zeta)*dx(q) = f(q, u) + D*dx(dx(q)) + sigma*q*eta
dt(u) + u*dx(u)        = g(q, u)

f(q, u) = q*(r + u - U0 - (r + delta)*(q-1)**2)
g(q, u) = epsilon1*(U0 -u) + epsilon2*(Ubar - u)*q

first, I do this without noise - sigma = 0

-Pratik Aghor
"""
#######################################
Nt = 10000;
nsave = Nt/100;
dt = 1.0e-2;

Lx = 200.0
Nx = 512

x_basis = de.Fourier('x', Nx, interval=(-Lx/2, Lx/2), dealias=3/2)
domain = de.Domain([x_basis], np.float64)

x = domain.grid(0)
#######################################

#######################################
make_dedalus_movie(Nt, nsave, Lx, Nx, [-1, 2], varname="u", name="u_movie", show=0)
make_dedalus_movie(Nt, nsave, Lx, Nx, [-10, 10], varname="q", name="q_movie", show=0)

#######################################
