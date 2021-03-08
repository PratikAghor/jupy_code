from numpy import pi, cos, arange, ones
import numpy as np
import params
from params import *
import cheb
from cheb import *
###########################################
'''
Inviscid SRI linear stability
author: Pratik Aghor
Ref: Viscous and inviscid strato-rotational instabilities, Robins et. al., JFM 2020.
'''
###########################################
###########################################
A00 = np.zeros((Nr+1, Nr+1), dtype=complex); A01 = A00; A02 = A00; A03 = A00; A04 = A00;

A10 = A00; A11 = A00; A12 = A00; A13 = A00; A14 = A00;

A20 = A00; A21 = A00; A22 = A00; A23 = A00; A24 = A00;

A30 = A00; A31 = A00; A32 = A00; A33 = A00; A34 = A00;

A40 = A00; A41 = A00; A42 = A00; A43 = A00; A44 = A00;
###########################################

B00 = A00; B01 = A00; B02 = A00; B03 = A00; B04 = A00;

B10 = A00; B11 = A00; B12 = A00; B13 = A00; B14 = A00;

B20 = A00; B21 = A00; B22 = A00; B23 = A00; B24 = A00;

B30 = A00; B31 = A00; B32 = A00; B33 = A00; B34 = A00;

B40 = A00; B41 = A00; B42 = A00; B43 = A00; B44 = A00;

###########################################

###########################################
"""
Find and sort the eigvals for each k
First define the easy rhs matrix B that does not depend on k
"""
B00 = D0r; B11 = D0r; B22 = D0r; B33=D0r;


"""
implement BCs u_r = 0.0 at r = r_i and r_o
We implement it as ur = -999 ur, because we want the 0 eigenvalue to
correspond to the critical k, m, Re; and not because of boundary conditions
"""
B00[0, 0:Nr+1] = -999.0; B00[Nr, 0:Nr+1] = -999.0;

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

    A00 = -1j*m*Omega*Id
    A01 = 2.0*Omega*Id - 2.0*1j*m*prefactor*oneByr*oneByr*Id
    A04 = D1r

    A10 = 2.0*1j*m*prefactor*oneByr*oneByr*Id - Z*Id
    A11 = -1j*m*Omega*Id
    A14 = -1j*m*oneByr*Id

    A22 = -1j*m*Omega*Id
    A23 = -Id
    A24 = -1j*k*Id

    A32 = oneByFrsq*Id
    A33 = -1j*m*Omega*Id

    A40 = D1r + oneByr*Id
    A41 = 1j*m*oneByr*Id
    A42 = 1j*k*Id

    # implement BCs u_r = 0.0 at r = r_i and r_o, no penetration
    A00[0, 0] = 1.0; A00[0, 1:Nr+1] = 0.0;
    A00[Nr, 0:Nr] = 0.0; A00[Nr, Nr] = 1.0;

    # stack them together
    A = np.block([ \
    [A00, A01, A02, A03, A04],
    [A10, A11, A12, A13, A14],
    [A20, A21, A22, A23, A24],
    [A30, A31, A32, A33, A34],
    [A40, A41, A42, A43, A44]
    ])
###########################################
