from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from numpy import pi, cosh, exp, round, zeros, identity, arange, real, cos, sin, multiply, transpose
from numpy.fft import fft,ifft, fftfreq
from matplotlib.pyplot import figure

import numpy as np  # Import numpy
import math
from numpy.linalg import inv
###########################################
'''
Code to test fft

The code is written on the transformed co-ordinates xc = [0, 2\pi]
Physical coordinates [0, 1]

- Pratik Aghor
'''
###########################################
# Grid:
# c stands for computational domain
xc_min = 0;
xc_max = 2*pi;

Lx = xc_max - xc_min;
Nx = 8;
dx = 1.0 / (Nx-1);

xc = (1.0*Lx/Nx)*arange(0, Nx);

samplingFreq = dx;
# u = sin(xc); # Initial condition
u = 2.0*sin(xc) + sin(3.0*xc); # Initial condition

v = fft(u)/len(u); # normalized fft

# v = v[range(int(len(u)/2))] # Exclude sampling frequency

print "v = fft(u) = \n", v
