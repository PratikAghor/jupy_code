from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from numpy import pi, cosh, exp, round, zeros, identity, arange, real, cos, sin, multiply, transpose
from numpy.fft import fft,ifft, fft2, ifft2
from matplotlib.pyplot import figure

import numpy as np  # Import numpy
import math
from numpy.linalg import inv
import h5py # for saving data

import heat_params
from heat_params import *

###########################################
'''
Code to integrate 2d heat equation (u_t = (\nabla^{2})^{2}u) with periodic
boundaries using Fourier spectral method.

For time integration, the second-order (two step) Adams-Moulton method for the linear term
and the second-order (two-step) Adams-Bashforth method for the nonlinear term are used (Henceforth called AM2AB2).

The code is written on the transformed co-ordinates (xc, yc) = [0, 2\pi] x [0, 2\pi]
Physical coordinates [0, Lx] x [0, Ly]

Throughout this code: index i is reserved for the y-direction, index j is reserved for the x-direction
and index n is reserved for time-marching

Nxc = Nx - 1; Nyc = Ny -1; are computational grid-points. We do not include the last point
because fft2 assumes periodic bc's in both directions. We pad the matrix later while saving
the solution. All computations are performed on (Nyc x Nxc) grid.

- Pratik Aghor
'''
###########################################
# Initial condition u0
u0 = zeros((Nyc, Nxc))
for i in range(0, Nyc):
    for j in range(0, Nxc):
        u0[i, j] = exp( -(x[j] - Lx/2.0)**2 - (y[i] - Ly/2.0)**2) # sin(xc[j]/2);  
# print "u(x, y) = sin(xc) \n"

u1 = u0;

# alpha_x = complex wavenumber in x direction
# alpha_y = complex wavenumber in y direction

alpha_x = 1j*kx;
alpha_y = 1j*ky;
"""
# Dx = d/d(xc) operator in Fourier space
# Dy = d/d(yc) operator in Fourier space
"""
Dx = a1*(alpha_x);
Dy = a2*(alpha_y);

"""
# L = linear operator in Fourier space
# v is the 2d fft of u

In fact, L should have been a 3d matrix with dimensions [Nyc, Nxc, Nxc],
in the sense, that for each y, we should invert have an Nxc x Nxc matrix.

However, in this particular case, the (Nxc x Nxc) matrix we would get is diagonal,
as the linear operator only consists of (Dx[j]*Dx[j] + Dy[i]*Dy[i]).

The A matrix from AM2AB2 algorithm (see documentation), is in fact a matrix of dimensions
[Nyc, Nxc, Nxc]. For each i value (for a fixed y), we have a diagonal matrix with constants.

In order to invert this matrix, we need only do inversion of the diagonal entries.

Here I write L[i, j] as a 2d matrix meant as an alias for the original 3d matrix.

Each i^(th) row of L[i, j] corresponds to the diagonal of the (Nxc x Nxc) matrix.

We must NOT treat L[i, j] as a regular 2d matrix, but think of each row of L[i, j]
to be a diagonal matrix. The A-matrix = (I - 0.5*dt*L) is then merely -

A[i, j] = (1.0 - 0.5*dt*L[i, j]).

The A-matrix again, is an alias for the original 3d matrix whose rows must be thought
of as diagonals of Nx x Nx matrices at each fiexed i.

Hence, in order to find Ainv, we need only invert each entry of A element-wise.

Similarly, we can create the matrix B[i, j] = (1.0 + 0.5*dt*L[i, j]).

Again, B[i, j] is a 2d stand-in for the original 3d matrix of size [Nyc, Nxc, Nxc].

By construction, the matmul(B, u) for example would become element-wise multiplication
multiply(B, u)

We must take care of this during time marching
"""
L = zeros((Nyc, Nxc), dtype=complex);
Ainv = zeros((Nyc, Nxc), dtype=complex); B = zeros((Nyc, Nxc), dtype=complex);

for i in range(0, Nyc):
    for j in range(0, Nxc):
        L[i, j]     = kappa*((Dx[j]*Dx[j] + Dy[i]*Dy[i]));
        Ainv[i, j]  = 1.0/(1.0 - 0.5*dt*L[i, j])
        B[i, j]     = (1.0 + 0.5*dt*L[i, j])

v0 = fft2(u0);
v1 = v0;

u1_padded = zeros((Ny, Nx)) # for plotting
#######################################
# time marching!
for n in range (0, Nt+1):

    u0 = real(ifft2( v0 ))
    u1 = real(ifft2( v1 ))

    v2 = multiply(Ainv, multiply(B, v1))

    # write u_padded with the last term = u[0] on Nx grid. for plotting purposes
    # padding in y
    u1_padded[0:Ny-1, 0:Nx-1] = u1[0:Ny-1, :];
    u1_padded[Ny-1, 0:Nx-1] = u1[0, :];

    # padding in x
    u1_padded[0:Ny-1, 0:Nx-1] = u1[:, 0:Nx-1];
    u1_padded[0:Ny-1, Nx-1] = u1[:, 0];

    # padding the corner
    u1_padded[Nx-1, Nx-1] = u1_padded[Nx-1, 0];

    # save data after every nsave steps
    if ((n % nsave) == 0):
        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('heat_data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('u(t)', data=u1_padded)
        hf.close()

    # update
    v1 = v2;
    v0 = v1;
