#######################################
# Compact finite difference schemes with specral-like resolution, Lele, J. Comp. Phys. 1992
# A compact FDM on staggered grid for Navier-Stokes flows, Zhang et. al., Int. J. Numer. Meth. Fluids. 2006;
# Author: Pratik Aghor
import numpy as np
import interpolate
from interpolate import *

import params
from params import *

#Import plotting functions:
import matplotlib as mpl
import matplotlib.pyplot as plt
#######################################
def lap_2d(f, x, y):
    """
    discrete 2d Laplacian f_xx + f_yy
    x, y, f are assumed to have periodic bcs
    f[i, j] = ith row, jth col
    """
    Ny, Nx = np.shape(f)

    f_2x = np.zeros((Ny, Nx))
    f_2y = np.zeros((Ny, Nx))

    for i in range(0, Ny):
        f_2x[i, :] = interpolate_colloc(f[i, :], x, A2x_colloc_pbc, B2x_colloc_pbc)

    for j in range(0, Nx):
        f_2y[:, j] = interpolate_colloc(f[:, j], y, A2y_colloc_pbc, B2y_colloc_pbc)

    lap_f = f_2x + f_2y

    return lap_f

# Schnakenberg nonlinearities
def f(u, v):
    return a - u + (u*u)*v

def g(u, v):
    return b - (u*u)*v
#######################################
if __name__ == "__main__":
    #This block will be evaluated if this script is called as the main routine
    #and will be ignored if this file is imported from another script.

    method = "AB2"
    # initialize
    u0 = np.zeros((ny, nx)); v0 = np.zeros((ny, nx))
    for i in range(0, ny):
        for j in range(0, nx):
            u0[i, j] = np.cos(x[j])*np.sin(y[i])
            v0[i, j] = np.cos(y[i])*np.sin(x[j])

    u1 = u0; v1 = v0;

    u_full = np.zeros((ny+1, nx+1)); v_full = np.zeros((ny+1, nx+1)) # for plotting


    # time marching using AB2
    for nt in range(0, Nt+1):

        lap_u0 = lap_2d(u0, x, y)
        lap_u1 = lap_2d(u1, x, y)

        lap_v0 = lap_2d(v0, x, y)
        lap_v1 = lap_2d(v1, x, y)

        # time-march u1, v1 equations
        if (method = "AB2"):
            u2 = u1 + (1.5*dt)*(lap_u1 + f(u1, v1)) - (0.5*dt)*(lap_u0 + f(u0, v0))
            v2 = v1 + (1.5*dt)*(d*lap_v1 + g(u1, v1)) - (0.5*dt)*(d*lap_v0 + g(u0, v0))

        if (nt % nsave == 0):
            # use pbc to fill out full u and v, save for post-processing
            u_full[0:ny, 0:nx] = u1
            v_full[0:ny, 0:nx] = v1

            u_full[:, nx] = u_full[:, 0]
            u_full[ny, :] = u_full[0, :]
            v_full[:, nx] = v_full[:, 0]
            v_full[ny, :] = v_full[0, :]

            np.savetxt('data/u'+str(nt)+'.asc', u_full)
            np.savetxt('data/v'+str(nt)+'.asc', v_full)

            print("step - ", nt, " done")

        # update for the next time-step
        u0 = u1
        u1 = u2
        v0 = v1
        v1 = v2


    #######################################
