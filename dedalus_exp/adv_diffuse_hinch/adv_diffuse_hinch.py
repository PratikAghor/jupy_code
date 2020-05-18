"""
1D advection-diffusion

dt(f) + dx(omega(x)*f) = \epsilon fxx
periodic BC on x in [0, 2\pi]
omega(x) is a known functon which is positive and periodic on x.

Ref: Perturbation Methods by Hinch Sec. # 7.3

This script should be ran serially (because it is 1D)

There should already be a folder named data

usage: in dedalus mode:
python3 adv_diffuse_hinch.py

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
# Bases and domain
Lx = 2.0*np.pi
Nx = 128
x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)
domain = de.Domain([x_basis], np.float64)

# Problem

# Problem
problem = de.IVP(domain, variables=['u', 'ux', 'uxx', 'omega', 'omegax'])

problem.parameters['epsilon'] = 0.01

#problem.add_equation("dt(u) + cf1*dx(ux) + cf2*dx(uxxx) = -cf3*ux*uxx")
problem.add_equation("dt(u) - epsilon*uxx = - omega*ux - u*omegax")
problem.add_equation("ux - dx(u) = 0")
problem.add_equation("uxx - dx(ux) = 0")
problem.add_equation("omega = 2.0 + sin(x)")
problem.add_equation("omegax - dx(omega) = 0")

# Build solver
solver = problem.build_solver(de.timesteppers.SBDF2)
solver.stop_wall_time = np.inf
solver.stop_iteration = 10000
nsave = solver.stop_iteration/500

# Initial conditions
x = domain.grid(0)
u = solver.state['u']
ux = solver.state['ux']
uxx = solver.state['uxx']
#################################################
# IC
n = 20; sigma = 2
u['g'] = np.sin(x)# np.log(1 + np.cosh(n)**2/np.cosh(n*(x-0.5*Lx))**2) / (2*n)
u.differentiate(0, out=ux)
ux.differentiate(0, out=uxx)

#################################################
# Store data for final plot
clean_data(); # clean files
u.set_scales(1)
write_data(0, np.copy(u['g']))

# Main loop
dt = 2e-3

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
