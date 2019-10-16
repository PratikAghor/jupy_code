from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import LineCollection
from numpy import pi, cosh, exp, round, zeros, identity, arange, real, cos, sin, multiply, transpose
from numpy.fft import fft,ifft
from matplotlib.pyplot import figure

import numpy as np  # Import numpy
import math
from numpy.linalg import inv
###########################################
'''
Code to integrate 1d wave u_t + c_{0}u_x=0 equation with periodic boundaries using Fourier spectral method. Crank-Nicholson CN was used. 
- Pratik Aghor
'''
###########################################
# Grid and required matrices:
Lx = 2*pi;
nsave = 5;
Nt = 400;
c_0 = 1.0; #wave-speed of u_t + c_{0}u_x=0
theta = 0.5; #weight to the current time-step value for the linear operator. theta = 0 => implicit, theta = 1 =>explicit, theta = 0.5 => Crank-Nicholson 
Nx = 256; dt = 0.01; x = (1.0*Lx/Nx)*arange(0, Nx);
#print "x=", x
u = exp(-(x-(Lx/2.0))*(x-(Lx/2.0))); #Initial condition

kx = zeros(Nx); kx[0:Nx/2] = arange(0,Nx/2); kx[Nx/2+1:] = arange(-Nx/2+1,0,1);
#print "kx=", kx
alpha = 2.0*pi*kx/Lx;              # real wavenumbers:    exp(alpha*x)    
#print "alpha=", alpha

D = (1j*alpha); # D = d/dx operator in Fourier space
#print "D =", '\n', D;

L = -c_0*D;                            
#linear operator in Fourier space, diagonal matrix!
#print "L=", L

A_inv = (np.ones(Nx) - (1.0-theta)*dt*L)**(-1);
#print "A_inv = ", '\n', A_inv
B     = np.ones(Nx) + theta*dt*L;

v  = fft(u);                  # transform u to spectral
#time marching!
for n in range (1, Nt):
    v = A_inv*(B*v);
   

    #Import plotting functions:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    fig = plt.figure()  # Create a figure instance
    ax = fig.gca()  
    ax.plot(x, real(ifft(v)))  # Plot the solution
    ax.set_xlabel('x')  # Set x label
    ax.set_ylabel('u')  # Set y label
    if n % nsave == 0:
        plt.savefig("wave_images/waveCN_"+str(n/nsave)+".png")
    #plt.show()  # Show the figure
    

