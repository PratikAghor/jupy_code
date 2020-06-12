"""
time marching for the Barkley excitable-bistable model
q = turbulence intensity
u = max velocity in the pipe direction (high for laminar, small for transition)
(as for puffs and slugs, the profile gets blunted)

dt(q) + (u-zeta)*dx(q) = f(q, u) + D*dx(dx(q)) + sigma*q*eta
dt(u) + u*dx(u)        = g(q, u)

f(q, u) = q*(r + u - U0 - (r + delta)*(q-1)**2)
g(q, u) = epsilon1*(U0 -u) + epsilon2*(Ubar - u)*q

first, I do this without noise - sigma = 0

Refs:
Theoretical perspective on the route to turbulence in a pipe - JFM Perspectives

This script should be ran serially (because it is 1D)

There should already be a folder named data

usage: in dedalus mode:
python3 barkley_turb_model.py

Pratik Aghor
"""
import h5py
import numpy as np
import time
import matplotlib.pyplot as plt
import os

from dedalus import public as de
from dedalus.extras.plot_tools import quad_mesh, pad_limits

import logging
logger = logging.getLogger(__name__)
#################################################
def clean_data():
    path = "data"
    fileList = os.listdir(path)
    for fileName in fileList:
        os.remove(path+"/"+fileName)
#################################################
def write_data(n, f, varname="u"):
    # save data after every nsave steps
    n_str = str(n)
    filename = str(varname) + n_str
    hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
    hf.create_dataset(str(varname), data=f)
    hf.close()
#################################################
#################################################
def write_ss_for_auto(x, u, ux, uxx, uxxx, fname='shper_1'): #
    # write the steady state in .dat file for auto 07p to read
    # assuming u, ux, uxx, uxxx have been scaled to account for dealiasing
    filename = str(fname) + '.dat'
    shperOut = np.column_stack((x, u))
    shperOut = np.column_stack((shperOut, ux))
    shperOut = np.column_stack((shperOut, uxx))
    shperOut = np.column_stack((shperOut, uxxx))
    np.savetxt(filename, shperOut)
#################################################
# Bases and domain
Lx = 200.0
Nx = 512
qc = 0.5

x_basis = de.Fourier('x', Nx, interval=(-Lx/2, Lx/2), dealias=3/2)
domain = de.Domain([x_basis], np.float64)

# Problem
problem = de.IVP(domain, variables=['q', 'qx', 'qxx', 'u', 'ux', 'uxx'])

problem.parameters['r'] = 0.6
problem.parameters['zeta'] = 0.8
problem.parameters['D'] = 0.5
problem.parameters['delta'] = 0.1
problem.parameters['epsilon1'] = 0.1
problem.parameters['epsilon2'] = 0.2
problem.parameters['U0'] = 1.0
problem.parameters['UBar'] = 0.5

problem.substitutions['f(q, u)'] = "q*(r + u - U0 - (r + delta)*(q-1)**2)"
problem.substitutions['g(q, u)'] = "epsilon1*(U0 - u) + epsilon2*(UBar - u)*q"

problem.add_equation("dt(q) - D*dx(qx) - zeta*qx =  f(q, u) - u*qx" )
problem.add_equation("dt(u) = - u*ux + g(q, u) ")
# auxillary equations
problem.add_equation("qx - dx(q) = 0")
problem.add_equation("qxx - dx(qx) = 0")
problem.add_equation("ux - dx(u) = 0")
problem.add_equation("uxx - dx(ux) = 0")

# Build solver
solver = problem.build_solver(de.timesteppers.SBDF2)
solver.stop_wall_time = np.inf
solver.stop_iteration = 10000
nsave = solver.stop_iteration/100

# Initial conditions
x = domain.grid(0)
q = solver.state['q']
qx = solver.state['qx']
qxx = solver.state['qxx']

u = solver.state['u']
ux = solver.state['ux']
uxx = solver.state['uxx']
#################################################
# IC
# localized initial condition (-> shper_1.dat)
# u['g'] = 10.0*(np.exp(-(x)**2/5.0**2))*np.cos(qc*(x))
# localized initial condition (-> shper_2.dat)
q['g'] = 1.0
u['g'] = -10.0*(np.exp( -(x)**2/5.0**2))*np.cos(3.0*(x))
q.differentiate(0, out=qx)
qx.differentiate(0, out=qxx)

u.differentiate(0, out=ux)
ux.differentiate(0, out=uxx)

#################################################
# Store data for final plot
clean_data(); # clean files
q.set_scales(1);u.set_scales(1);
write_data(0, np.copy(q['g']), varname="q")
write_data(0, np.copy(u['g']))

# Main loop
dt = 1.0e-2

try:

    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        solver.step(dt)
        if solver.iteration % nsave == 0:
            q.set_scales(1);u.set_scales(1);
            # c_list.append(np.copy(c['g']))
            # t_list.append(solver.sim_time)
            write_data(solver.iteration, np.copy(q['g']), varname="q")
            write_data(solver.iteration, np.copy(u['g']))

        if solver.iteration % 100 == 0:
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))

#################################################
# write steady state for auto-07p
# set scales to account for dealiasing
# ux.set_scales(1); uxx.set_scales(1); uxxx.set_scales(1);
#
# write_ss_for_auto(x, np.copy(u['g']), np.copy(ux['g']),\
#  np.copy(uxx['g']), np.copy(uxxx['g']), fname='shper_2')
#################################################
