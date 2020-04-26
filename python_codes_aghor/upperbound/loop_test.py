from matplotlib.collections import LineCollection
from numpy import pi, arange, cos, sin, zeros, identity, matmul, multiply
from numpy.fft import fft,ifft
from matplotlib.pyplot import figure

import numpy as np  # Import numpy
import math
from numpy.linalg import inv, solve
import h5py # for saving data

###########################################
'''
loop test
n represnts time-, i represnts z-, j represnts x-coordinate

Author: Pratik Aghor
'''
###########################################
###########################################
nx = 8
ny = 4
for i in range(0, ny):
    for j in range(0, nx):
        print "in loop # 1 \n"

for i in range(0, ny):
    for j in range(0, nx):
        print "in loop # 2 \n"


###########################################
