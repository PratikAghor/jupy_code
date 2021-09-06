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
    f is assumed to have periodic bc
    f is also assumed to have one ghost node in each direction
    """
    Ny, Nx = np.shape(f)

    f_2x = np.zeros((Ny, Nx))
    f_2y = np.zeros((Ny, Nx))

    for i in range(0, Ny):
        f_2x[i, :] = D2x @ f[i, :]

    for j in range(0, Nx):
        f_2y[:, j] = D2y @ f[:, j]

    lap_f = f_2x + f_2y

    return lap_f

# Schnakenberg nonlinearities
def f(u, v):
    return a - u + (u*u)*v

def g(u, v):
    return b - (u*u)*v

def apply_pbc(u, nx, ny):
    u[:, nx+1]  = u[:, 1]
    u[:, 0]     = u[:, nx]
    u[ny+1, :]  = u[1, :]
    u[0, :]     = u[ny, :]
#######################################
if __name__ == "__main__":
    #This block will be evaluated if this script is called as the main routine
    #and will be ignored if this file is imported from another script.

    method = "AB2"
    # initialize
    u0 = np.zeros((ny+2, nx+2)); v0 = np.zeros((ny+2, nx+2))
    for i in range(1, ny+1):
        for j in range(1, nx+1):
            u0[i, j] = np.cos(x[j-1])*np.sin(y[i-1])
            v0[i, j] = np.cos(y[i-1])*np.sin(x[j-1])

    # pbc
    apply_pbc(u0, nx, ny)
    apply_pbc(v0, nx, ny)

    u1 = u0; v1 = v0;

    u_full = np.zeros((ny+1, nx+1)); v_full = np.zeros((ny+1, nx+1)) # for plotting


    # time marching using AB2
    for nt in range(0, Nt+1):

        lap_u0 = lap_2d(u0, x, y)
        lap_u1 = lap_2d(u1, x, y)

        lap_v0 = lap_2d(v0, x, y)
        lap_v1 = lap_2d(v1, x, y)

        # time-march u1, v1 equations
        if (method == "AB2"):
            u2 = u1 + (1.5*dt)*(lap_u1 + f(u1, v1)) - (0.5*dt)*(lap_u0 + f(u0, v0))
            v2 = v1 + (1.5*dt)*(d*lap_v1 + g(u1, v1)) - (0.5*dt)*(d*lap_v0 + g(u0, v0))

        if (nt % nsave == 0):
            # use pbc to fill out full u and v, save for post-processing
            u_full = u1[1:, 1:]
            v_full = v1[1:, 1:]

            np.savetxt('data/u'+str(nt)+'.asc', u_full)
            np.savetxt('data/v'+str(nt)+'.asc', v_full)

            print("step - ", nt, " done")

        # update for the next time-step
        u0 = u1
        u1 = u2
        v0 = v1
        v1 = v2
        apply_pbc(u0, nx, ny)
        apply_pbc(v0, nx, ny)


    #######################################
