"""
time marching for 1d Swift-Hohenberg equation (sh1d)

dt(u) = ru - [(dx**2 + qc**2)**2] u + v*u**2-g*u**3

Refs:
1. Localized states in the generalized Swift-Hohenberg equation
by John Burke and Edgar Knobloch

2. Münsteranian Torturials on Nonlinear Science
by Uwe Thiele, Oliver Kamps, Svetlana Gurevich

This script should be ran serially (because it is 1D)

There should already be a folder named data

usage: in dedalus mode:
python3 sh1d.py

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
problem = de.IVP(domain, variables=['u', 'ux', 'uxx', 'uxxx'])

problem.parameters['q_c'] = qc
problem.parameters['v'] = 0.41
problem.parameters['g'] = 1.0
problem.parameters['r'] = -0.01245
problem.parameters['uhom'] = 0.0

#problem.add_equation("dt(u) + cf1*dx(ux) + cf2*dx(uxxx) = -cf3*ux*uxx")
problem.add_equation("dt(u) - r*u + (dx(uxxx) + 2.0*(q_c**2)*uxx +(q_c**4)*u)  = v*u**2-g*u**3")
problem.add_equation("ux - dx(u) = 0")
problem.add_equation("uxx - dx(ux) = 0")
problem.add_equation("uxxx - dx(uxx) = 0")

# Build solver
solver = problem.build_solver(de.timesteppers.SBDF2)
solver.stop_wall_time = np.inf
solver.stop_iteration = 20000
nsave = solver.stop_iteration/500

# Initial conditions
x = domain.grid(0)
u = solver.state['u']
ux = solver.state['ux']
uxx = solver.state['uxx']
uxxx = solver.state['uxxx']
#################################################
# IC
# localized initial condition (-> shper_1.dat)
# u['g'] = 10.0*(np.exp(-(x)**2/5.0**2))*np.cos(qc*(x))
# localized initial condition (-> shper_2.dat)
u['g'] = -10.0*(np.exp( -(x)**2/5.0**2))*np.cos(qc*(x))
u.differentiate(0, out=ux)
ux.differentiate(0, out=uxx)
uxx.differentiate(0, out=uxxx)

#################################################
# Store data for final plot
clean_data(); # clean files
u.set_scales(1)
write_data(0, np.copy(u['g']))

# Main loop
dt = 1.0e-2

try:

    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        solver.step(dt)
        if solver.iteration % nsave == 0:
            u.set_scales(1)
            # c_list.append(np.copy(c['g']))
            # t_list.append(solver.sim_time)
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
ux.set_scales(1); uxx.set_scales(1); uxxx.set_scales(1);

write_ss_for_auto(x, np.copy(u['g']), np.copy(ux['g']),\
 np.copy(uxx['g']), np.copy(uxxx['g']), fname='shper_2')
#################################################
