from numpy import pi, cos, arange, ones

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

import sys
sys.path.append('/home/aghor/UNH/independent/jupy_codes_local/dedalus_exp/sri/sri_lin_stab')

# for saving figs in this folder
import os
my_path = os.path.abspath(__file__)

from functions import *
from functions.params import *
from functions.cheb import *

###########################################
'''
SRI linear stability
author: Pratik Aghor
'''
###########################################
A00 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A01 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
A02 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A03 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A04 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

A10 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A11 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
A12 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A13 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A14 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

A20 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A21 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
A22 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A23 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A24 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

A30 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A31 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
A32 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A33 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A34 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

A40 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A41 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
A42 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A43 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A44 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);
###########################################

B00 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B01 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
B02 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B03 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B04 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

B10 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B11 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
B12 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B13 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B14 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

B20 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B21 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
B22 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B23 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B24 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

B30 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B31 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
B32 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B33 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B34 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

B40 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B41 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); \
B42 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B43 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex); B44 = A00 = np.zeros((Nr+1, Nr+1), dtype=complex);

###########################################
# scalar Laplacian
def Lap_s(k):
    Lap_s = D2r + np.multiply(oneByr, D1r) \
- np.multiply( (m**2+1.0)*(np.multiply(oneByr, oneByr)) + k**2 , D0r)
    return Lap_s
###########################################
"""
Find and sort the eigvals for each k
First define the easy rhs matrix B that does not depend on k
"""
B00 = D0r; B11 = D0r; B22 = D0r; B33=D0r;
# stack them together
B = np.block([ \
[B00, B01, B02, B03, B04],
[B10, B11, B12, B13, B14],
[B20, B21, B22, B23, B24],
[B30, B31, B32, B33, B34],
[B40, B41, B42, B43, B44]
])

for k_counter in range(0, Nk):
    k = kvec[k_counter]
    Lap = Lap_s(k)

    A00 = 1j*m*Omega*D0r + (eta/(1.0-eta))*(1.0/Re)*Lap
    A01 = np.multiply( (2.0*Omega - 2.0*1j*m*np.multiply(oneByr, oneByr)), D0r)
    A04 = - D1_g2gl_c

    A10 = np.multiply((- Z + 2.0*1j*m*np.multiply(oneByr, oneByr)), D0r)
    A11 = 1j*m*Omega*D0r + (eta/(1.0-eta))*(1.0/Re)*Lap
    A14 = np.multiply(-1j*m*oneByr, D0_g2gl)

    A22 = 1j*m*Omega*D0r + (eta/(1.0-eta))*(1.0/Re)*Lap
    A23 = -D0r
    A24 = (-1j*k)*D0_g2gl

    A32 = oneByFrsq*D0r
    A33 = -1j*m*Omega*D0r

    A40 = D1r + np.multiply(oneByr, D0r)
    A41 = np.multiply(1j*m*oneByr, D0r)
    A42 = np.multiply(1j*k, D0r)

    # implement BCs u_r = u_theta = u_z = 0.0 at r = r_i and r_o
    A00[0, 0] = 1.0; A00[0, 1:Nr+1] = 0.0; A01[0, 0:Nr+1] = 0.0; A02[0, 0:Nr+1] = 0.0;\
    A03[0, 0:Nr+1] = 0.0; A04[0, 0:Nr+1] = 0.0;

    A00[Nr, 0:Nr] = 0.0; A00[Nr, Nr] = 1.0; A01[Nr, 0:Nr+1] = 0.0; A02[Nr, 0:Nr+1] = 0.0;\
    A03[Nr, 0:Nr+1] = 0.0; A04[Nr, 0:Nr+1] = 0.0;

    A10[0, 0:Nr+1] = 0.0; A11[0, 0] = 1.0; A11[0, 1:Nr+1] = 0.0; A12[0, 0:Nr+1] = 0.0;\
    A13[0, 0:Nr+1] = 0.0; A14[0, 0:Nr+1] = 0.0;

    A10[Nr, 0:Nr+1] = 0.0; A11[Nr, 0:Nr] = 0.0; A11[Nr, Nr] = 1.0; A12[Nr, 0:Nr+1] = 0.0;\
    A13[Nr, 0:Nr+1] = 0.0; A14[Nr, 0:Nr+1] = 0.0;

    A20[0, 0:Nr+1] = 0.0; A21[0, 0:Nr+1] = 0.0; A22[0, 0] = 1.0; A22[0, 1:Nr+1] = 0.0;\
    A23[0, 0:Nr+1] = 0.0; A24[0, 0:Nr+1] = 0.0;

    A20[Nr, 0:Nr+1] = 0.0; A21[Nr, 0:Nr+1] = 0.0; A22[Nr, 0:Nr] = 0.0; A22[Nr, Nr] = 1.0;\
    A23[Nr, 0:Nr+1] = 0.0; A24[Nr, 0:Nr+1] = 0.0;

    # stack them together
    A = np.block([ \
    [A00, A01, A02, A03, A04],
    [A10, A11, A12, A13, A14],
    [A20, A21, A22, A23, A24],
    [A30, A31, A32, A33, A34],
    [A40, A41, A42, A43, A44]
    ])

    """
    # TODO:
    (1) Solve the eigenvalue problem Ax = sigma_c Bx
    (2) Perhaps it is beneficial to define \rho at the Gauss-Grid
    """
