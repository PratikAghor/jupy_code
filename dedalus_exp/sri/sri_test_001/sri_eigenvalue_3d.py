"""
Usage:
python3 sri_eigenvalue_3d.py [--re=<reynolds> --eta=<eta> --G=<grashof> --epsilon=<epsilon> --Pr=<prandtl> --m=<initial_m> --k=<k> --mu=<mu> --ar=<Gamma> --alpha_z=<alpha_z>]

Options:
  --re=<inner_reynolds>  Reynolds number for simulation [default: 5182]
  --eta=<eta>      Eta - ratio of R1/R2 [default: 0.6]
  --G=<grashof>    Grashof number [default: 50]
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

nr = 32

args=docopt(__doc__)
Re1=float(args['--re'])
eta=np.float(args['--eta'])
G=float(args['--G'])
epsilon = float(args['--epsilon'])
Pr = float(args['--Pr'])
m = int(args['--m'])
k = int(args['--k'])
mu = np.float(args['--mu'])
Gamma = float(args['--ar'])
alpha_z = float(args['--alpha_z'])

"""
delta = R2 - R1
mu = Omega2/Omega1
eta = R1/R2

scale [L] = delta
scale [T] = delta/(R1 Omega1)
scale [V] = R1 Omega1

Trying to implement The Boussinesq approximation in rapidly rotating flows (2013, JFM)
"""
#derived parameters
R1 = eta/(1. - eta)
R2 = 1./(1-eta)
midpoint = R1 + (R2-R1) / 2
if alpha_z:
    Lz = 2*np.pi/alpha_z
else:
    Lz = Gamma

logger.info("Re:{:.3e}, eta:{:.4e}, G:{:.3e}, epsilon:{:.3e}, Pr:{:.4e}, mu:{:.4e}".format(Re1,eta,G,epsilon,Pr,mu))

logger.info("Lz set to {:.6e}".format(Lz))

variables = ['u','ur','v','vr','w','wr','p'] #'T','Tr',

#domain
r_basis = de.Chebyshev('r', nr, interval=[R1, R2])

bases = [r_basis]
domain = de.Domain(bases)

#problem
problem = de.EVP(domain, eigenvalue='lambda', variables=variables)

#params into equations
problem.parameters['eta']=eta
problem.parameters['mu']=mu
problem.parameters['Lz']=Lz
problem.parameters['R1']=R1
problem.parameters['R2']=R2
problem.parameters['Re_i']=Re1
problem.parameters['Re_o']=mu*Re1
problem.parameters['Pr']=Pr
problem.parameters['epsilon']=epsilon
problem.parameters['pi']=np.pi
problem.parameters['kz'] = 2*np.pi/Lz * k
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

problem.substitutions['A'] = '(Re_o - eta*Re_i)/(1.0 + eta)'
problem.substitutions['B'] = 'eta*(Re_i - eta*Re_o)/((1.0-eta)*(1.0-eta**2))'
problem.substitutions['C'] = '-(4.0*np.log(eta) + (1.0-eta**2)*(3.0-eta**2) )\
/(16.0*(1.0-eta**2)*(( (1.0+eta**2)*np.log(eta) + (1.0-eta**2) )) )'

problem.substitutions['v0'] = 'A*r + B/r'       #background profile? forcing instead of forcing the boundaries
problem.substitutions['dv0dr'] = 'A - B/(r*r)'  #d/dr of background forcing

problem.substitutions['w0'] = 'G*(C*(r**2 - r_i**2) + (C*(r_o**2 - r_i**2) + 0.25*(r_o**2 - r**2))*np.log(r/r_i)/np.log(eta) )'
problem.substitutions['dw0dr'] = 'G*(C*(2*r) + 0.25*(-r_o**2/r - 2*r*np.log(r/r_i) - r)/np.log(eta) )'

problem.substitutions['T0'] = '0.5 + np.log(r/r_i)/np.log(eta)'
problem.substitutions['dT0dr'] = '1.0/(r*np.log(eta))'
problem.substitutions['dp0dr(Tb, vb)'] = '1.0 - epsilon*Tb*vb*vb/r'  #d/dr of background pressure
problem.substitutions['dp0dz'] = '4.0*C + 0.5 -  1.0/np.log(eta)'  #d/dz of background pressure

problem.substitutions['dtheta(f)'] = '1j*m*f'
problem.substitutions['dz(f)'] = '1j*kz*f'
problem.substitutions['dt(f)'] = 'lambda*f'

# assume pre-multiplication by r*r
problem.substitutions['Lap_s(f, f_r)'] = "r*r*dr(f_r) + r*f_r + dtheta(dtheta(f)) + r*r*dz(dz(f))"
problem.substitutions['Lap_r'] = "Lap_s(u, ur) - u - 2*dtheta(v)"
problem.substitutions['Lap_t'] = "Lap_s(v, vr) - v + 2*dtheta(u)"
problem.substitutions['Lap_t_v0'] = "r*r*dr(dv0dr) + r*dv0dr - v0"
problem.substitutions['Lap_z'] = "Lap_s(w, wr)"
problem.substitutions['Lap_z_w0'] = "r*r*dr(dw0dr) + r*dw0dr"


# momentum equations
problem.add_equation("r*r*dt(u) - Lap_r + (1 - epsilon*T0)*(-2*r*v0*v + r*v0*dtheta(u) + r*r*w0*dz(u))+ r*r*dr(p) = -r*r*dp0dr(T0, v0)")
problem.add_equation("r*r*dt(v) - Lap_t + (1 - epsilon*T0)*(r*r*dv0dr*u + r*v0*u + r*v0*dtheta(v) + r*r*w0*dz(v) )+ r*dtheta(p)  = Lap_t_v0")
problem.add_equation("r*r*dt(w) - Lap_z + r*r*dz(p) + (1 - epsilon*T0)*(r*v0*dtheta(w) + r*r*w0*dz(w))= Lap_z_w0")

# energy equation
# problem.add_equation("r*r*dt(T) - (1.0/Pr)*Lap_s(T, Tr) + (r*v0*dtheta(T) + r*r*w0*dz(T) + u*dT0dr= (1.0/Pr)*Lap_s(T0, dT0dr)")

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
