"""
Dedalus script for 2D reaction-diffusion equations.

This script uses a Chebyshev basis in the y-direction with Neumann boundary
conditions and Fourier in x-direction with periodic BC's

Apparently, dedalus cannot use Chebyshev in two directions!

Ref: Subcritical Turing bifurcation and the morphogenesis of localized patterns,
Vı́ctor Breña-Medina and Alan Champneys, PRE, 2014

This script can be ran serially or in parallel, and uses the built-in analysis
framework to save data snapshots in HDF5 files.  The `merge_procs` command can
be used to merge distributed analysis sets from parallel runs, and the
`plot_slices.py` script can be used to plot the snapshots.

To run, merge, and plot using 4 processes, for instance, you could use:
    $ mpiexec -n 4 python3 rayleigh_benard.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs snapshots
    $ mpiexec -n 4 python3 plot_slices.py snapshots/*.h5

This script can restart the simulation from the last save of the original
output to extend the integration.  This requires that the output files from
the original simulation are merged, and the last is symlinked or copied to
`restart.h5`.

To run the original example and the restart, you could use:
    $ mpiexec -n 4 python3 react_diffuse_sys_1.py
    $ mpiexec -n 4 python3 -m dedalus merge_procs snapshots
    $ ln -s snapshots/snapshots_s2.h5 restart.h5
    $ mpiexec -n 4 python3 react_diffuse_sys_1.py

The simulations should take a few process-minutes to run.

"""

import numpy as np
from mpi4py import MPI
import time
import pathlib

from dedalus import public as de
from dedalus.extras import flow_tools

import logging
logger = logging.getLogger(__name__)


# Parameters
Lx, Ly = (500, 500)

# Create bases and domain
x_basis = de.Fourier('x', 256, interval=(-Lx/2, Lx/2), dealias=3/2)
y_basis = de.Chebyshev('y', 256, interval=(-Ly/2, Ly/2))
domain = de.Domain([x_basis, y_basis], grid_dtype=np.float64)

# 2D uniaxial Swift-Hohenberg (SH)
problem = de.IVP(domain, variables=['u', 'ux', 'uxx', 'uy', 'uyy', 'v', 'vx', 'vxx', 'vy', 'vyy'])

problem.parameters['r'] = 1.0e-4
problem.parameters['b'] = 1.0e-4
problem.parameters['c'] = 0.1
problem.parameters['h'] = 1.0e-2
problem.parameters['k2'] = 4.0e-3
problem.parameters['D1'] = 0.13
problem.parameters['D2'] = 10.0

problem.add_equation("dt(u) - D1*(uxx + uyy) + (c + r)*u - h*v     =  k2*u*u*v")
problem.add_equation("dt(v) - D2*(vxx + vyy) - c*u       + h*v     = -k2*u*u*v + b")

problem.add_equation("ux    - dx(u)     = 0")
problem.add_equation("uxx   - dx(ux)    = 0")
problem.add_equation("uy    - dy(u)     = 0")
problem.add_equation("uyy   - dy(uy)    = 0")
problem.add_equation("vx    - dx(v)     = 0")
problem.add_equation("vxx   - dx(vx)    = 0")
problem.add_equation("vy    - dy(v)     = 0")
problem.add_equation("vyy   - dy(vy)    = 0")

logger.info("Neumann BC: fixed flux in y direction, periodic in x")

problem.add_bc( "left(uy) = 0")
problem.add_bc("right(uy) = 0")

problem.add_bc( "left(vy) = 0")
problem.add_bc("right(vy) = 0")


# Build solver
solver = problem.build_solver(de.timesteppers.SBDF3)
logger.info('Solver built')

# Initial conditions or restart
if not pathlib.Path('restart.h5').exists():

    # Initial conditions
    x = domain.grid(0, scales=problem.domain.dealias)
    y = domain.grid(1)

    u = solver.state['u']; ux = solver.state['ux']; uy = solver.state['uy'];
    uxx = solver.state['uxx']; uyy = solver.state['uyy'];
    v = solver.state['v']; vx = solver.state['vx']; vy = solver.state['vy'];
    vxx = solver.state['vxx']; vyy = solver.state['vyy'];

    u.differentiate('x',  out=ux)
    ux.differentiate('x',  out=uxx)
    u.differentiate('y',  out=uy)
    uy.differentiate('y',  out=uyy)

    v.differentiate('x',  out=vx)
    vx.differentiate('x',  out=vxx)
    v.differentiate('y',  out=vy)
    vy.differentiate('y',  out=vyy)

    # IC
    u.set_scales(domain.dealias, keep_data=False)
    u.set_scales(domain.dealias, keep_data=False)

    # u['g'] = 5.0*(np.exp( -(x)**2 - (y)**2 ))
    u['g'] = np.sin(7*x/Lx)*np.sin(7*y/Ly)
    v['g'] = 5.0*np.sin(3*x/Lx)*np.sin(3*y/Ly)
    # Timestepping and output
    dt = 0.5
    stop_sim_time = 500
    fh_mode = 'overwrite'

else:
    # Restart
    write, last_dt = solver.load_state('restart.h5', -1)

    # Timestepping and output
    dt = last_dt
    stop_sim_time = 200
    fh_mode = 'append'

# Integration parameters
solver.stop_sim_time = stop_sim_time

# Analysis
snapshots = solver.evaluator.add_file_handler('snapshots', sim_dt=0.25, max_writes=50, mode=fh_mode)
snapshots.add_system(solver.state)

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.proceed:
        dt = 0.5
        dt = solver.step(dt)
        if (solver.iteration-1) % 10 == 0:
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
