#######################################
# Computation of Strongly Compressible Rotating Flows, Harada, J. Comp. Phys., 1980;
# A compact FDM on staggered grid for Navier-Stokes flows, Zhang et. al., Int. J. Numer. Meth. Fluids. 2006;
# Author: Pratik Aghor
#######################################
import numpy as np
from numpy.linalg import *

from functions.cheb import *
from functions.params import *
from functions.enforce_bc import *
#######################################

#Import plotting functions:
import matplotlib as mpl
import matplotlib.pyplot as plt
#######################################
# check if import params is ok
# print "Ekman # = ", E
#######################################
def get_fu_rhs(u, v, w, rho, T, p, rhoBar, rhoBaru, rhoBarv, rhoBarw, rhoBarT):
	"""
	rhs for r-momentum equation
	"""
	urhoBaru = np.multiply(u, rhoBaru)
	wrhoBaru = np.multiply(w, rhoBaru)

	for i in range(0, Nz):
		urhoBaru_1r[i, :] = np.matmul(D1r, urhoBaru[i, :])

	for j in range(0, Nr):
		urhoBaru_1z[:, j] = np.matmul(D1z, wrhoBaru[:, j])
	#---------------------------------
	for i in range(0, Nzc):
		for j in range(0, Nrf):
			fu_rhs[i, j] = -(urhoBaru_1r_rf[i, j] + wrhoBaru_1z_rf[i, j]) \
			+ rhoBar[i, j] * rc[j] + rhoBarv[i, j] * (2.0 + vByr[i, j]) \
			- oneBygammaMsq * rc[j] * p_1r_rf[i, j]\
			+ E * rc[j] * ( (4.0/3.0) * (u_2r[i, j] + uByr_1r[i, j])\
			+ u_2z[i, j] + (1.0/3.0) * w_rz[i, j] );
	#---------------------------------
	# if abs value less than cutoff, set to zero
	# fu_rhs[np.abs(fu_rhs) < cutoff] = 0
	#---------------------------------
	return fu_rhs
#######################################
#######################################
def get_fv_rhs(u, v, w, rho, T, p, rhoBar, rhoBaru, rhoBarv, rhoBarw, rhoBarT):
	"""
	rhs for phi-momentum equation
	"""

	#---------------------------------
	# if abs value less than cutoff, set to zero
	# fv_rhs[np.abs(fv_rhs) < cutoff] = 0
	#---------------------------------
	return fv_rhs
#######################################
#######################################
def get_fw_rhs(u, v, w, rho, T, p, rhoBar, rhoBaru, rhoBarv, rhoBarw, rhoBarT):
	"""
	rhs for z-momentum equation
	"""
	#---------------------------------
	# if abs value less than cutoff, set to zero
	# fw_rhs[np.abs(fw_rhs) < cutoff] = 0
	#---------------------------------
	return fw_rhs
#######################################
#######################################
def get_frho_rhs(rhoBaru, rhoBarw):
	# fields required for get_frho_rhs

	#---------------------------------
	for i in range(0, Nzc):
		rhoBaru_1r_c[i, :] = interpolate_f2c(rhoBaru[i, :], rf, A1r_f2c, B1r_f2c);

	for j in range(0, Nrc):
		rhoBarw_1z_c[:, j] = interpolate_f2c(rhoBarw[:, j], zf, A1z_f2c, B1z_f2c);
	#---------------------------------
	frho_rhs = -(rhoBaru_1r_c + rhoBarw_1z_c);
	#---------------------------------
	# if abs value less than cutoff, set to zero
	# frho_rhs[np.abs(frho_rhs) < cutoff] = 0
	#---------------------------------
	return frho_rhs
#######################################
#######################################
def get_fT_rhs(u, v, w, rho, T, p, rhoBar, rhoBaru, rhoBarv, rhoBarw, rhoBarT):
	#---------------------------------

	#---------------------------------
	#---------------------------------
	# Q = \nabla . u
	#---------------------------------
	for i in range(0, Nzc):
		for j in range(0, Nrc):
			Q[i, j] = (1.0/rc[j])*(ru_1r_c[i, j]) + w_1z_c[i, j];
	#---------------------------------
	#---------------------------------
	# build Phi
	#---------------------------------
	#---------------------------------
	# DEBUG:
	#---------------------------------
	for i in range(0, Nzc):
		for j in range(0, Nrc):
			Phi[i, j] = mu*( 2.0*u_1r_c[i, j]*u_1r_c[i, j] \
			+ ((2.0/(rc[j]*rc[j])) * u_c[i, j]*u_c[i, j])\
			+ 2.0*w_1z_c[i, j]*w_1z_c[i, j] \
			+ (v_1r_c[i, j] - (v_c[i, j]/rc[j]))*(v_1r_c[i, j] - (v_c[i, j]/rc[j]))\
			+(u_1z_c[i, j] + w_1r_c[i, j])*(u_1z_c[i, j] + w_1r_c[i, j])\
			+ v_1z_c[i, j]*v_1z_c[i, j]\
			- (2.0/3.0)*Q[i, j]*Q[i, j] );

	#---------------------------------
	# print("PHI = \n ", Phi)

	#---------------------------------
	# build \nabla^2 (T)
	#---------------------------------
	# use T_rf, T_zf to encorporate bc's effect on derivatives
	#---------------------------------
	for i in range(0, Nzc):
		for j in range(0, Nrc):
			fT_rhs[i, j] = -(rhoBaru_1r_c[i, j] + rhoBarw_1z_c[i, j])*T[i, j]\
			- (gamma-1.0)*rhoBarT[i, j]*Q[i, j]\
			+ gammaEbyPr*rc[j]*(T_2r_c[i, j] + (1.0/rc[j])*T_1r_c[i, j] + T_2z_c[i, j] )\
			+ (gamma-1.0)*gamma*M*M*E*rc[j]*Phi[i, j];

	#---------------------------------
	# if abs value less than cutoff, set to zero
	# fT_rhs[np.abs(fT_rhs) < cutoff] = 0
	#---------------------------------
	return fT_rhs
#######################################
