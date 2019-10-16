###################################################
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from numpy import exp
import sympy as sp
from numpy.random import rand, random
###################################################

#parameters
v_m = 0.05;#; #max car velocity
rho_m = 0.05; #max car density

R = 0.1; # R<rho_m IC = R exp(-x^2/L^2)
L = 1.0;

tb_Bar = np.sqrt(np.exp(1))/(2.*np.sqrt(2)) * (L/R) * (rho_m/v_m);
print "time of incipient breaking is ", tb_Bar
#the grid
x_min = -4.
x_max = 4.
N = 100 # #of grid points
tmax = tb_Bar + 1;


dt  = 0.01;
Nt  = int (tmax/dt);
#print "Nt =", Nt
#x0 = np.linspace(x_min, x_max, N)  # the 1d domain
x = np.linspace(x_min, x_max, N)  # the 1d domain
#print "x =", x

rho0 = np.zeros(N);
c0 = np.zeros(N);
rho = np.zeros((Nt + 1, N))


t   = 0.0; # set t = 0 for every node before beginning

for i in range(Nt):
    t = i *  dt;
    for j in range(N):
        rho0[j] = R * exp (-(x[j]/L)**2);

        def f(x0):
            return x[j] -  v_m * (1.0 - (2.*R/rho_m) * exp(-(x0/L)**2))* t - x0;
            # as c0 = v_m * (1.0 - (2.*R/rho_m) * exp(-(x0/L)**2));
            #arbitrary initial guess
        x0 = fsolve(f, x[j] - v_m *  (1.0 - (2.*R/rho_m) * exp(-(x[j]/L)**2)) * t )
        #print "t = ", t, " x = ", x[j], "x0 =", x0
        rho[i, j] = R * exp (-(x0/L)**2);


fig = plt.figure()  # Create a figure instance
ax = fig.gca()  # Get current axes
ax.plot(x, rho0)  # Plot the solution
for i in range(Nt+1):
    if i % 200 == 0:
        ax.plot(x, rho[i])
ax.set_xlabel('x')  # Set x label
ax.set_ylabel(r'$\rho_{0}$')  # Set y label
plt.savefig("traffic_flow_method_of_characteristics_solution.png")
plt.show()  # Show the figure
