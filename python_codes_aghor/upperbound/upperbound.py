from matplotlib.collections import LineCollection
from numpy import pi, arange, cos, sin, zeros, identity, matmul, multiply
from numpy.fft import fft,ifft
from matplotlib.pyplot import figure

import numpy as np  # Import numpy
import math
from numpy.linalg import inv, solve
import h5py # for saving data

import params
from params import *

import cheb
from cheb import *

import qrsolve
from qrsolve import *
###########################################
'''
Code to integrate Euler-Lagrange equations for porous medium convection
with Dirichlet boundaries using Chebyeshev spectral method.
For time integration, the second-order Adams-Moulton (Crank-Nicholson) method for the linear term
and the second-order (two-step) Adams-Bashforth method for the nonlinear term are used (Henceforth called AM2AB2).

The code is written on physical coordinates z = [0, 1], x = [0, 2]

Cheb derivatives are taken on the Gauss-Lobatto grid zc = [-1, 1]

Ref: Computational approaches to aspect-ratio-dependent upper bounds
and heat flux in porous medium convection - Wen et. al., Physics Letters A

n represnts time-, i represnts z-, j represnts x-coordinate

Author: Pratik Aghor
'''
###########################################
###########################################
# declare variables and IC
theta0 = zeros((Nz+1, Nx)); tau0 = zeros(Nz+1);
W0 = zeros((Nz+1, Nx)); gamma0 = zeros((Nz+1, Nx));

theta1 = zeros((Nz+1, Nx)); tau1 = zeros(Nz+1);
W1 = zeros((Nz+1, Nx)); gamma1 = zeros((Nz+1, Nx));

theta2 = zeros((Nz+1, Nx)); tau2 = zeros(Nz+1);
W2 = zeros((Nz+1, Nx)); gamma2 = zeros((Nz+1, Nx));

I = identity(Nz+1);
Wj_thetaj_0 = zeros((Nz+1, Nx)); Wj_thetaj_1 = zeros((Nz+1, Nx));

# IC
for i in range(0, Nz+1):
    tau0[i] = 1.0 - z[i];
    for j in range(0, Nx):
        theta0[i, j] = sin((j+1)*pi*z[i]);

tau1 = tau0; theta1 = theta0;

Ainv = zeros((Nz+1, Nz+1, Nx));

for j in range(0, Nx):
    k = kvec[j];
    A0 = ((a**2)*D2 - (k**2)*I); # the (a**2) comes due to cheb grid- coordinate transform from z to zc

    # Dirichlet BC
    A0[0, 0] = 1.0; A0[0, 1:Nz+1] = 0.0;
    A0[Nz, 0:Nz] = 0.0; A0[Nz, Nz] = 1.0;

    Ainv[:, :, j] = inv(A0);

for j in range(0, Nx):
    k = kvec[j];
    W0[:, j] = (-ra * k * k) * matmul(Ainv[:, :, j], theta0[:, j]);
    W1[:, j] = W0[:, j];

    gamma0[:, j] = matmul(Ainv[:, :, j], multiply( (-ra * a * matmul(D1, tau0)), theta0[:, j]) );

    gamma1[:, j] = gamma0[:, j];
###########################################
# declare rhs matrices, vecs, etc
# A = zeros((Nz+1, Nz+1)); overwrite the above A for each k as required
b = zeros(Nz+1);
###########################################
# time marching
for n in range(0, Nt):
    for j in range(0, Nx):
        k = kvec[j];

        # step # 1: update theta
        Ltheta = 2.0 * ((a**2)*D2 - (k**2)*identity(Nz+1));

        Ntheta0 = multiply( -(a * matmul(D1, tau0)), W0[:, j] ) - k * k * gamma0[:, j];
        Ntheta1 = multiply( -(a * matmul(D1, tau1)), W1[:, j] ) - k * k * gamma1[:, j];

        # solve A * theta2 = b problem to update theta
        A = I - (0.5 * dt)*Ltheta;
        b = matmul((I + (0.5 * dt)*Ltheta) , theta1[:, j]) + 0.5 * dt*(3.0* Ntheta1 - Ntheta0);

        # Dirichlet BC
        A[0, 0] = 1.0; A[0, 1:Nz+1] = 0.0;
        A[Nz, 0:Nz] = 0.0; A[Nz, Nz] = 1.0;
        b[0] = 0.0; b[Nz] = 0.0;

        theta2[:, j] = solve(A, b);
        # print "size(theta2) = ", np.size(theta2)

        # step # 2: update W
        b = (-ra * k * k)*theta2[:, j]
        # Dirichlet BC
        b[0] = 0.0; b[Nz] = 0.0;

        W2[:, j] = matmul(Ainv[:, :, j], b); #qrsolve(A0, b); #

        # step # 3: solve for tau
        Wj_thetaj_0[:, j] = multiply(W1[:, j], theta1[:, j]);
        Wj_thetaj_1[:, j] = multiply(W2[:, j], theta2[:, j]);
        # end the 1st j-loop here

    # sum Wj_thetaj_0 column-wise
    sum_Wj_thetaj_0 = zeros(Nz+1); sum_Wj_thetaj_1 = zeros(Nz+1);
    for i in range(0, Nz+1):
        for j in range(0, Nx):
            sum_Wj_thetaj_0[i] = sum_Wj_thetaj_0[i] + Wj_thetaj_0[i, j];
            sum_Wj_thetaj_1[i] = sum_Wj_thetaj_1[i] + Wj_thetaj_1[i, j];

    Ltau = (a**2) * D2;
    A = (I - 0.5 * dt * Ltau);
    b = matmul( (I + 0.5 * dt * Ltau), tau1 ) \
    - 0.125*dt*(3.0*a*matmul(D1, sum_Wj_thetaj_1) - a*matmul(D1, sum_Wj_thetaj_0));

    # Dirichlet BC
    A[0, 0] = 1.0; A[0, 1:Nz+1] = 0.0;
    A[Nz, 0:Nz] = 0.0; A[Nz, Nz] = 1.0;
    b[0] = 0.0; b[Nz] = 1.0;
    # print "\n A for tau = \n", A
    # print "\n b for tau = \n", b

    tau2 = solve(A, b);

    # step # 4: solve for gamma
    for j in range(0, Nx):
        # solve A * gamma2 = b problem to update gamma2
        b = multiply( (-ra * a * matmul(D1, tau2)), theta2[:, j]);

        # Dirichlet BC
        b[0] = 0.0; b[Nz] = 0.0;

        gamma2[:, j] = matmul(Ainv[:, :, j], b); # qrsolve(A0, b); #
        # end second j-loop here


    # save data after every nsave steps
    if ((n % nsave) == 0):
        n_str = str(n)
        filename = "tau" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('tau(t)', data=tau1)
        hf.close()

        filename = "W" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('W(t)', data=W1)
        hf.close()

        filename = "gamma" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('gamma(t)', data=gamma1)
        hf.close()

        filename = "theta" + n_str
        hf = h5py.File('data/' + filename +'.h5', 'w') # open file in h5 format
        hf.create_dataset('theta(t)', data=theta1)
        hf.close()

    # update all variables for time-marching
    theta0 = theta1;
    theta1 = theta2;

    W0 = W1;
    W1 = W2;

    tau0 = tau1;
    tau1 = tau2;

    gamma0 = gamma1;
    gamma1 = gamma2;

