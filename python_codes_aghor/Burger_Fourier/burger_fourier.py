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
Code to integrate 1d Burger's equation (u_t + u*u_x = nu*u_xx) with periodic boundaries using Fourier spectral method. Crank-Nicholson, Adams-Bashforth second order accurate scheme (CNAB2) was used for time stepping. 

- Pratik Aghor
'''
###########################################
# Grid and required matrices:
Lx = 4*pi;
nsave = 50;
Nt = 5000;
nu = 0.1;
theta = 0.5; #weight to the current time-step value for the linear operator. theta = 0 => implicit, theta = 1 =>explicit, theta = 0.5 => Crank-Nicholson 
Nx = 256; dt = 0.002; x = (1.0*Lx/Nx)*arange(0, Nx);
#print "x=", x
u = exp(-(x-(Lx/2.0))*(x-(Lx/2.0))); #Initial condition
#print "u=", u;
#print "u*u=", u*u;

kx = zeros(Nx); kx[0:Nx/2] = arange(0,Nx/2); kx[Nx/2+1:] = arange(-Nx/2+1,0,1);
#print "kx=", kx
alpha = 2.0*pi*kx/Lx;              # real wavenumbers:    exp(alpha*x)    
#print "alpha=", alpha

D = (1j*alpha);# D = d/dx operator in Fourier space
#print "D =", '\n', D;

L = nu*D*D;                                     # linear operator nu*D^2 in Fourier space, diagonal matrix!
#print "L=", L
G = -0.5*D;                                     # -1/2 D operator in Fourier space        

Nsave = Nt/nsave + 1;                           # number of saved time steps, including t=0
t = zeros(Nsave);
t[0:Nsave] = arange(0, Nsave)*(dt*nsave);       # t timesteps


U = zeros((Nsave, Nx));                         # matrix of (ti, xj) values.
U[0,:] = u;                                     # assign initial condition to U
counter = 2 ;                                   # counter for saved data


# convenience variables
dt2  = 0.5*dt;
dt32 = 1.5*dt;
#print dt32

#A     = (np.ones(Nx) - dt2*L)
#print "A = ", '\n', A
A_inv = (np.ones(Nx) - (1.0-theta)*dt*L)**(-1);
#print "A_inv = ", '\n', A_inv
B     = np.ones(Nx) + theta*dt*L;

N1  = G*fft(u*u);                               #-u u_x (spectral), notation N1 = N^n = N(u(n dt))
N0 = N1;                                        #notation N0 = N^{n-1} = N(u((n-1) dt))
v  = fft(u);                                    # transform u to spectral

#time marching!
for n in range (1, Nt):
    N0 = N1;
    N1 = G*fft((real(ifft( v )))*(real(ifft( v )))); # compute N0 = -u u_x
    

    v = A_inv*((B*v) + dt32*N1 - dt2*N0);
    
    if n % nsave == 0:
        U[counter, :]  = real(ifft(v))
        counter = counter + 1

    #Import plotting functions:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    fig = plt.figure()  # Create a figure instance
    ax = fig.gca()  
    ax.plot(x, real(ifft(v)))  # Plot the solution
    ax.set_xlabel('x')  # Set x label
    ax.set_ylabel('u')  # Set y label
    if n % nsave == 0:
        plt.savefig("burger_images/u_"+str(n/nsave)+".png")
    #plt.show()  # Show the figure
    
#print "U =", '\n', U;
#Import plotting functions:
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#from matplotlib import cm

#plt.pcolor(x,t,U)
#plt.xlabel("$x$")
#plt.ylabel("$t$")
#plt.xlim((0, Lx))
#cbar = plt.colorbar()
#plt.tight_layout()
#plt.show()

#from matplotlib import cm

# create the figure
#fig = plt.figure()
#ax = fig.add_subplot(111)
#im = ax.imshow(U, cmap=cm.coolwarm)
#cbar = fig.colorbar(im)
#plt.show()
