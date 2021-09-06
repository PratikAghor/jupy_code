#######################################
# Author: Pratik Aghor
import numpy as np
from numpy.linalg import *
#######################################
#######################################
def get_d1x_mat(n, dh):
	"""
	2d mat for getting f_x
	one ghost node in each dir assumed
	for example
	f_1x[i, :] = D*f_i if f_i is the ith row
	len(f_i) = n, including ghosts
	"""
	D = np.zeros((n, n))

	for i in range(1, n-1):
		D[i, i-1] = -1
		D[i, i+1] = 1
	D = (1/(2*dh))*D

	return D
#######################################
#######################################
def get_d2x_mat(n, dh):
	"""
	2d mat for getting f_xx
	one ghost node in each dir assumed
	for example
	f_2x[i, :] = D2*f_i if f_i is the ith row
	len(f_i) = n, including ghosts
	"""
	D2 = np.zeros((n, n))

	for i in range(1, n-1):
		D2[i, i-1] = 1
		D2[i, i]   = -2
		D2[i, i+1] = 1
	D2 = (1/(dh*dh))*D2

	return D2
#######################################
#######################################
def get_d1y_mat(n, dh):
	"""
	2d mat for getting f_y
	one ghost node in each dir assumed
	for example
	f_1x[:, j] = D*f_j if f_j is the jth col
	len(f_j) = n, including ghosts

	NOTE: the index i and the co-ordinate y run opposite
	"""
	Dy = np.zeros((n, n))

	for i in range(1, n-1):
		Dy[i, i-1] = 1
		Dy[i, i+1] = -1
	Dy = (1/(2*dh))*Dy

	return Dy
#######################################
#######################################
def get_d2y_mat(n, dh):
	"""
	2d mat for getting f_y
	one ghost node in each dir assumed
	for example
	f_2x[:, j] = D2y*f_j if f_j is the jth col
	len(f_j) = n, including ghosts

	NOTE: the index i and the co-ordinate y run opposite
	"""
	D2y = np.zeros((n, n))

	for i in range(1, n-1):
		D2y[i, i+1] = 1
		D2y[i, i]   = -2
		D2y[i, i-1] = 1
	D2y = (1/(dh*dh))*D2y

	return D2y
#######################################
