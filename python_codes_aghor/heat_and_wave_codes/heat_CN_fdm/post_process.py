# post processing
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import params
from params import *
#######################################
"""
post-process heat_CN data
"""
#######################################

sigma_exact = -nu*(4.0*np.pi)*(4.0*np.pi);

# declare an array that stores the value of u at a point
log_value_at_a_point = np.zeros((Nt/nsave) + 1);
tsave = np.zeros((Nt/nsave) + 1);

fig = plt.figure(1)  # Create a figure instance
ax = fig.gca()

for n in range(0, Nt):
    # construct the exact solution
    u_exact = np.exp(sigma_exact*t[n])*np.sin(4.0*(np.pi)*x);
    # read data
    if ((n % nsave) == 0):

        n_str = str(n)
        filename = "u" + n_str
        hf = h5py.File('data/' + filename + '.h5', 'r')
        ut = hf.get('u(t)')
        u = np.array(ut)

        log_value_at_a_point[n/nsave] = math.log(u[10]/ (float)(u_0[10]));
        tsave[n/nsave] = t[n]

        ax.plot(x, u_exact, '-', linewidth=2.0, label = r'$u_{exact}$' )  # Plot the exact solution
        ax.plot(x, u, 'o', linewidth=2.0, fillstyle='none', markeredgewidth = 2, label = r'$u$')  # Plot the numerical solution
        ax.set_xlabel('x')  # Set x label

#ax.legend()
plt.savefig('ut.png')
#######################################
# calculate numerical sigma = dissipation rate
sigma = (log_value_at_a_point[6] - log_value_at_a_point[5])/(float)(tsave[6] - tsave[5]);

sigma_str = str(sigma)
sigma_exact_str = str(sigma_exact)

print "calculated dissipation rate from code = ",'\t', sigma
print "exact dissipation rate = ",'\t', sigma_exact

# plot log (u_at_a_point) vs t
fig2 = plt.figure(2) # Create a figure instance
ax2 = fig2.gca()
ax2.plot(tsave, log_value_at_a_point, '--', linewidth=2.0)
ax2.set_ylabel(r'$log(u[10]/u_{0}[10])$')  # Set y label
ax2.set_xlabel('t')  # Set x label

plt.text(0.7, 0.9,r'$\sigma_n$ = ' + sigma_str,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)

plt.text(0.7, 0.8,r'$\sigma_{exact}$ = ' + sigma_exact_str,
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax.transAxes)
plt.savefig('sigma_compare.png')
#######################################
