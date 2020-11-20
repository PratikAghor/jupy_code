"""
Usage:
python3 sri_eigenvalue_3d.py [--re=<reynolds> --eta=<eta> --G=<grashof> --epsilon=<epsilon> --pr=<prandtl> --m=<initial_m> --k=<k> --mu=<mu> --ar=<Gamma> --alpha_z=<alpha_z>]

Options:
  --re=<inner_reynolds>  Reynolds number for simulation [default: 5182]
  --eta=<eta>      Eta - ratio of R1/R2 [default: 0.6]
  --Fr=<grashof>    Froude number [default: 1]
  --epsilon=<epsilon> Relative density [default: 0]
  --Pr<prandtl>     Prandtl number [default: 4.35]
  --m=<initial_m>  m mode [default: 6]
  --k=<k>          k mode [default: 1]
  --mu=<mu>        mu = Omega2/Omega1 [default: 0]
  --ar=<Gamma>     Aspect ratio (height/width) [default: 3]
  --alpha_z=<alpha_z>  z wavenumber [default: 3.13]
"""

import numpy as np
from dedalus import public as de
import logging
from docopt import docopt
from eigentools import Eigenproblem
logger = logging.getLogger(__name__)
"""
Trying to implement Viscous and inviscid strato-rotational instability by Robins et. al. (2020, JFM)
"""

nr = 32
#args=docopt(__doc__)

# geometric parameters
# Gamma = 10.0; # aspect ratio - in Hyun and Park = H/L
eta = 0.9; # eta = r_i/r_o, radius ratio
mu = 0.4; # Omega_o/Omega_i
Re1 = 100;
Fr = 1.0;
oneByFrsq = 1.0/(Fr**2)
alpha_z=3.13

# other parameters
m = 1
k = 4

#derived parameters
Re2 = mu*Re1
R1 = eta/(1. - eta)
R2 = 1./(1-eta)
midpoint = R1 + (R2-R1) / 2
Lz = 2*np.pi/alpha_z

logger.info("Re:{:.3e}, eta:{:.4e}, Fr:{:.3e}, mu:{:.4e}".format(Re1,eta,Fr,mu))

logger.info("Lz set to {:.6e}".format(Lz))

variables = ['u','ur','v','vr','w','wr','rho','p'] #

#domain
r_basis = de.Chebyshev('r', nr, interval=[R1, R2])

bases = [r_basis]
domain = de.Domain(bases)

#problem
problem = de.EVP(domain, eigenvalue='sigma_c', variables=variables)

#params into equations
problem.parameters['eta']=eta
problem.parameters['mu']=mu
problem.parameters['Lz']=Lz
problem.parameters['R1']=R1
problem.parameters['R2']=R2
problem.parameters['Re_i']=Re1
problem.parameters['Re_o']=Re2
problem.parameters['Fr']=Fr
# problem.parameters['Pr']=Pr
# problem.parameters['epsilon']=epsilon
problem.parameters['pi']=np.pi
problem.parameters['kz'] = (2*np.pi/Lz) * k
problem.parameters['m'] = m

#Substitutions

"""
this implements the cylindrical del operators.
NB: ASSUMES THE EQUATION SET IS PREMULTIPLIED BY A POWER OF r (SEE BELOW)!!!

Lap_s --> scalar laplacian
Lap_r --> r component of vector laplacian
Lap_t --> theta component of vector laplacian
Lap_z --> z component of vector laplacian

"""
problem.substitutions['r_i'] = 'eta/(1.0 - eta)'
problem.substitutions['r_o'] = '1.0/(1.0 - eta)'

problem.substitutions['A'] = '(mu - eta*eta)/(1.0 - eta*eta)'
problem.substitutions['B'] = '(eta**2)*(1.0 - mu)/((1.0+eta)*(1.0-eta)**3)'
problem.substitutions['Z'] = '2.0*A'

problem.substitutions['Omega'] = 'A + B/(r*r)'       # background profile? forcing instead of forcing the boundaries

problem.substitutions['dtheta(f)'] = '1j*m*f'
problem.substitutions['dz(f)'] = '1j*kz*f'
problem.substitutions['dt(f)'] = 'sigma_c*f'

# assume pre-multiplication by r*r
problem.substitutions['Lap_s(f, f_r)'] = "r*r*dr(f_r) + r*f_r + dtheta(dtheta(f)) + r*r*dz(dz(f))"
problem.substitutions['Lap_r'] = "Lap_s(u, ur) - u - 2*dtheta(v)"
problem.substitutions['Lap_t'] = "Lap_s(v, vr) - v + 2*dtheta(u)"
problem.substitutions['Lap_z'] = "Lap_s(w, wr)"

# momentum equations
problem.add_equation("r*r*dt(u) - (eta/(1.0-eta))*(1.0/Re_i)*Lap_r + (-2*Omega*v*r*r + r*r*Omega*dtheta(u) )+ r*r*dr(p) = 0")
problem.add_equation("r*r*dt(v) - (eta/(1.0-eta))*(1.0/Re_i)*Lap_t + (r*r*Z*u + r*r*Omega*dtheta(v) )+ r*dtheta(p)  = 0")
problem.add_equation("r*r*dt(w) - (eta/(1.0-eta))*(1.0/Re_i)*Lap_z + r*r*dz(p) + r*r*Omega*dtheta(w) + r*r*rho = 0")

# energy equation (for rho - density perturbation)
problem.add_equation("r*r*dt(rho) + (r*r*Omega*dtheta(rho)) - (1.0/(Fr*Fr))*w = 0")

#continuity
problem.add_equation("r*ur + u + dtheta(v) + r*dz(w) = 0")

#Auxillilary equations
problem.add_equation("ur - dr(u) = 0")
problem.add_equation("vr - dr(v) = 0")
problem.add_equation("wr - dr(w) = 0")
# problem.add_equation("Tr - dr(T) = 0")

#Boundary Conditions
problem.add_bc("left(u) = 0")
problem.add_bc("right(u) = 0")
problem.add_bc("left(v) = 0")
problem.add_bc("right(v) = 0")
problem.add_bc("left(w) = 0")
problem.add_bc("right(w) = 0")
# problem.add_bc("left(T) = 0")
# problem.add_bc("right(T) = 0")


ep = Eigenproblem(problem)

growth, index, freq = ep.growth_rate({})

logger.info("Growth rate = {:16.15e}; frequency = {:16.15e}".format(growth, freq[0]))

#ep.spectrum(spectype='hires')
