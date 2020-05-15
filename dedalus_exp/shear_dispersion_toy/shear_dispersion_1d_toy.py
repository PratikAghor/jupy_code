"""
1D shear dispersion toy model

This script should be ran serially (because it is 1D)

there should already be a folder named data

usage: in dedalus mode:
python3 shear_dispersion_1d_toy.py

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
def write_data(n, f):
    # save data after every nsave steps
    n_str = str(n)
    filename = "c" + n_str
    hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
    hf.create_dataset('c(t)', data=f)
    hf.close()
#################################################
# Bases and domain
Lx = 2.0*np.pi
x_basis = de.Fourier('x', 128, interval=(-Lx, Lx), dealias=3/2)
domain = de.Domain([x_basis], np.float64)

# Problem
problem = de.IVP(domain, variables=['c', 'cx', 'cxx', 'cxxx', 'cxxxx', 'cxxxxx', 'cxxxxxx'])
problem.parameters['nu'] = 2e-4
problem.parameters['mu'] = 1e-4

problem.add_equation("dt(c) - nu*dx(cx) - dx(dx(cxxxx)) = 2.0*c*cxxxx + 4.0*cx*cxxx + 1.0 - mu*c")
problem.add_equation("cx - dx(c) = 0")
problem.add_equation("cxx - dx(cx) = 0")
problem.add_equation("cxxx - dx(cxx) = 0")
problem.add_equation("cxxxx - dx(cxxx) = 0")
problem.add_equation("cxxxxx - dx(cxxxx) = 0")
problem.add_equation("cxxxxxx - dx(cxxxxx) = 0")
# Build solver
solver = problem.build_solver(de.timesteppers.SBDF2)
solver.stop_wall_time = np.inf
solver.stop_iteration = 20000
nsave = solver.stop_iteration/1000

# Initial conditions
x = domain.grid(0)
c = solver.state['c']
cx = solver.state['cx']
cxx = solver.state['cxx']
cxxx = solver.state['cxxx']
cxxxx = solver.state['cxxxx']
cxxxxx = solver.state['cxxxxx']
cxxxxxx = solver.state['cxxxxxx']
#################################################
n = 20
c['g'] = np.log(1 + np.cosh(n)**2/np.cosh(n*x)**2) / (2*n)
c.differentiate(0, out=cx)
cx.differentiate(0, out=cxx)
cxx.differentiate(0, out=cxxx)
cxxx.differentiate(0, out=cxxxx)
cxxxx.differentiate(0, out=cxxxxx)
cxxxxx.differentiate(0, out=cxxxxxx)

#################################################
# Store data for final plot
clean_data(); # clean files
c.set_scales(1)
write_data(0, np.copy(c['g']))

# c_list = [np.copy(c['g'])]
# t_list = [solver.sim_time]

# Main loop
dt = 2e-3

try:

    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        solver.step(dt)
        if solver.iteration % nsave == 0:
            c.set_scales(1)
            # c_list.append(np.copy(c['g']))
            # t_list.append(solver.sim_time)
            write_data(solver.iteration, np.copy(c['g']))

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

# # Create space-time plot
# c_array = np.array(c_list)
# t_array = np.array(t_list)
# #xmesh, ymesh = quad_mesh(x=x, y=t_array)
# plt.figure()
# for n in range(0, len(t_array)):
#     if(n%10 == 0):
#         plt.plot(x, c_array[n, :])
#
# # plt.axis(pad_limits(xmesh, ymesh))
# # plt.colorbar()
# plt.xlabel('$x$')
# plt.ylabel('$u$')
# plt.title('toy_1d, (nu,mu)=(%g,%g)' %(problem.parameters['nu'], problem.parameters['mu']))
# plt.savefig('toy_1d.png')
