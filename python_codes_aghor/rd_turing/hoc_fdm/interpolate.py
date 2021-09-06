#######################################
# spin-up Hyun and Park JFM, 1992;
# A compact FDM on staggered grid for Navier-Stokes flows, Zhang et. al., Int. J. Numer. Meth. Fluids. 2006;
# Author: Pratik Aghor
import numpy as np
from numpy.linalg import *
import qrsolve
from qrsolve import *
#######################################
#######################################
# use a different function to create matrices
def create_matrices_interpolate_0_f2c(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 0th derivative (function value itself)
	from cell-face to cell-center.

    Inputs:
    x: grid (where f is defined), i.e., co-ordinates of cell-faces only

    Outputs:
    A, B to solve A varphi = B psi
	"""
	nface = len(x)

	n = nface - 1;

	# declare matrices
	A = np.zeros((n, n)); B = np.zeros((n, nface));

    # create matrices
	alpha = 3.0/10.0    # alpha for inside the domain
	a = (1.0/8.0)*(9.0 + 10.0*alpha)
	b = (1.0/8.0)*(6.0*alpha - 1.0)

	alpha_b = 1.0/4.0; # alpha for boundary scheme Nagarajan et al.
	a_b = 1.0/16.0*(5.0-alpha_b)
	b_b = 1.0/16.0*(9.0*alpha_b + 15.0)
	c_b = 1.0/16.0*(9.0*alpha_b - 5.0)
	d_b = 1.0/16.0*(1.0-alpha_b)
	# build A
	# Zhang et. al. boundary scheme - didn't work at the right boundary
	# A[0, 0] = 1.0; A[0, 1] = 1.0;

	# Nagarajan et al boundary scheme
	A[0, 0] = 1.0; A[0, 1] = alpha_b;


	for i in range(1, n-1):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;

	# Zhang et. al. boundary scheme - didn't work at the right boundary
	# A[n-1, n-2] = 1.0; A[n-1, n-1] = 1.0;

	# Nagarajan et al boundary scheme
	A[n-1, n-2] = alpha_b; A[n-1, n-1] = 1.0;

	# build B
	# Zhang et. al. boundary scheme - didn't work at the right boundary
	# B[0,0] = 1.0/4.0; B[0, 1] = 3.0/2.0; B[0, 2] = 1.0/4.0;

	# Nagarajan et al boundary scheme
	B[0,0] = a_b; B[0, 1] = b_b; B[0, 2] = c_b; B[0, 4] = d_b;
	B[1,0] = b/2.0; B[1, 1] = a/2.0; B[1, 2] = a/2.0; B[1, 3] = b/2.0;
	for i in range(2, n-2):
		B[i, i-1] = b/2.0;
		B[i, i]   = a/2.0;
		B[i, i+1] = a/2.0;
		B[i, i+2] = b/2.0;
	B[n-2,nface-4] = b/2.0; B[n-2, nface-3] = a/2.0; B[n-2, nface-2] = a/2.0; B[n-2, nface-1] = b/2.0;

	# Zhang et. al. boundary scheme - didn't work at the right boundary
	# B[n-1,n-3] = 1.0/4.0; B[n-1, n-2] = 3.0/2.0; B[n-1, n-1] = 1.0/4.0;

	# Nagarajan et al boundary scheme
	B[n-1,nface-4] = d_b; B[n-1,nface-3] = c_b; B[n-1, nface-2] = b_b; B[n-1, nface-1] = a_b;
	return A, B
#######################################

#######################################
# use a different function to create matrices
def create_matrices_interpolate_1_f2c(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 1st derivative
	from cell-face to cell-center.

    Inputs:
    x: grid (where f is defined), i.e., co-ordinates of cell-faces only

    Outputs:
    A, B to solve A varphi' = B psi
	"""
	nface = len(x)
	dx = (x[nface-1] - x[0])/(nface - 1)
	n = nface - 1;

	# declare matrices
	A = np.zeros((n, n)); B = np.zeros((n, nface));

    # create matrices
	alpha = 9.0/62.0    # alpha for inside the domain
	a = (3.0/8.0)*(3.0 - 2.0*alpha)
	b = (1.0/8.0)*(-1.0 + 22.0*alpha)
	# weights for discrete conservation
	w0 = 223.0/186.0
	w1 = 61.0/62.0

	alpha_b = 0.0; # alpha for boundary scheme Nagarajan et al. for discrete conservation
	a_b = (1.0/24.0)*(alpha_b - 23.0)
	b_b = (1.0/8.0)*(-9.0*alpha_b + 7.0)
	c_b = (1.0/8.0)*(9.0*alpha_b + 1.0)
	d_b = (-1.0/24.0)*(alpha_b + 1.0)

    # create matrices
	alpha_1 = 37.0/183.0    # alpha for inside the domain
	a1 = (3.0/8.0)*(3.0 - 2.0*alpha_1)
	b1 = (1.0/8.0)*(-1.0 + 22.0*alpha_1)

	# build A

	# Nagarajan et al boundary scheme with discrete conservation
	A[0, 0] = 1.0*w0; A[0, 1] = alpha_b*w0;
	A[1, 0] = alpha_1*w1; A[1, 1] = 1.0*w1; A[1, 2] = alpha_1*w1;

	for i in range(2, n-2):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;

	# Nagarajan et al boundary scheme with discrete conservation
	A[n-2, n-3] = alpha_1*w1; A[n-2, n-2] = 1.0*w1; A[n-2, n-1] = alpha_1*w1;
	A[n-1, n-2] = alpha_b*w0; A[n-1, n-1] = 1.0*w0;

	# build B
	# Nagarajan et al boundary scheme with discrete conservation
	B[0,0] = (a_b/dx)*w0; B[0, 1] = (b_b/dx)*w0; B[0, 2] = (c_b/dx)*w0; B[0, 3] = (d_b/dx)*w0;
	B[1,0] = (-b1/(3.0*dx))*w1; B[1, 1] = (-a1/dx)*w1; B[1, 2] = (a1/dx)*w1; B[1, 3] = (b1/(3.0*dx))*w1;
	for i in range(2, n-2):
		B[i, i-1] = -b/(3.0*dx);
		B[i, i]   = -a/dx;
		B[i, i+1] = a/dx;
		B[i, i+2] = b/(3.0*dx);
	# Nagarajan et al boundary scheme with discrete conservation
	B[n-2,nface-4] = (-b1/(3.0*dx))*w1; B[n-2, nface-3] = (-a1/dx)*w1; B[n-2, nface-2] = (a1/dx)*w1; B[n-2, nface-1] = (b1/(3.0*dx))*w1;
	B[n-1,nface-4] = (-d_b/dx)*w0; B[n-1,nface-3] = (-c_b/dx)*w0; B[n-1, nface-2] = (-b_b/dx)*w0; B[n-1, nface-1] = (-a_b/dx)*w0;
	return A, B
#######################################

#######################################
# use a different function to create matrices
def create_matrices_interpolate_2_f2c(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 2nd derivative
	from cell-face to cell-center.

    Inputs:
    x: grid (where f is defined), i.e., co-ordinates of cell-faces only

    Outputs:
    A, B to solve A varphi'' = B psi

	Inside domain: Zhang et al
	Boundary scheme: Zhang et al

	NOTE: viscous terms need to be evaluated on only collocated grids. This
	function might not be useful, but writing it for completeness.
	"""
	nface = len(x)
	dx = (x[nface-1] - x[0])/(nface - 1)
	n = nface - 1

	# declare matrices
	A = np.zeros((n, n)); B = np.zeros((n, nface));

    # create matrices
	alpha = 5.0/14.0    # alpha for inside the domain
	a = 6.0/7.0
	b = 6.0/7.0

	# boundary scheme Zhang et al.
	alpha_b = 1.0;
	a_b = 2.0
	b_b = -4.0
	c_b = 2.0

    # create matrices
	# build A

	# boundary scheme Zhang et al.
	A[0, 0] = 1.0; A[0, 1] = alpha_b;
	A[1, 0] = alpha; A[1, 1] = 1.0; A[1, 2] = alpha;

	for i in range(2, n-2):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;

	# boundary scheme Zhang et al.
	A[n-2, n-3] = alpha; A[n-2, n-2] = 1.0; A[n-2, n-1] = alpha;
	A[n-1, n-2] = alpha_b; A[n-1, n-1] = 1.0;

	# build B
	# boundary scheme Zhang et al.
	B[0, 0] = a_b; B[0, 1] = b_b; B[0, 2] = c_b;
	B[1, 0] = a; B[1, 1] = -2.0*a + b; B[1, 2] = a - 2.0*b; B[1, 3] = b;
	for i in range(2, n-2):
		B[i, i-1] = a;
		B[i, i]   = -2.0*a + b;
		B[i, i+1] = a - 2.0*b;
		B[i, i+2] = b;
	# boundary scheme Zhang et al.
	B[n-2,nface-4] = a; B[n-2, nface-3] = -2.0*a + b; B[n-2, nface-2] = a - 2.0*b; B[n-2, nface-1] = b;
	B[n-1,nface-3] = c_b; B[n-1, nface-2] = b_b; B[n-1, nface-1] = a_b;

	B = (1.0/(dx*dx))*B
	return A, B
#######################################
#######################################

def interpolate_f2c(f, x, A, B):

	"""
	Interpolate the 0th derivative (function value itself)
	from cell-face to cell-center.

	Inputs:
	f: discrete function
	x: grid (where f is defined), i.e., co-ordinates of cell-faces only
	A, B : obtained from create_matrices_interpolate_zero_f2c
	Outputs:
	f_at_center: interpolated f at cell-centers
	"""

	n = len(x)
	f_interpolated = np.zeros(n) # declare f_interpolated
	f_at_center	= np.zeros(n-1)

	rhs = np.matmul(B, f)

	f_at_center = qrsolve(A, rhs)

	# f_at_center = solve(A, rhs)
	return f_at_center
#######################################

# center-to-face (c2f) stuff
#######################################
#######################################
# use a different function to create matrices
def create_matrices_interpolate_0_c2f(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 0th derivative (function value itself)
	from cell-center to cell-face (c2f).

    Inputs:
    x: grid (where f is defined), i.e., co-ordinates of cell-centers only

    Outputs:
    A, B to solve A varphi = B psi
	"""
	nc = len(x)

	nf = nc + 1; # nface = nc + 1f

	# declare matrices
	A = np.zeros((nf, nf)); B = np.zeros((nf, nc));

    # create matrices
	alpha_b = 0.0; # alpha for boundary scheme Nagarajan et al.
	a_b = (1.0/8.0)*(3.0*alpha_b + 15.0)
	b_b = (1.0/4.0)*(3.0*alpha_b - 5.0)
	c_b = (1.0/8.0)*(3.0 - alpha_b)

	alpha_1 = 1.0/6.0    # alpha_1 for j = 1 and j = nc (second and second last to boundary)
	a_1 = (1.0/8.0)*(9.0 + 10.0*alpha_1)
	b_1 = (1.0/8.0)*(6.0*alpha_1 - 1.0)

	alpha = 3.0/10.0    # alpha for inside the domain
	a = (1.0/8.0)*(9.0 + 10.0*alpha)
	b = (1.0/8.0)*(6.0*alpha - 1.0)

	# build A
	# Nagarajan et al boundary scheme
	A[0, 0] = 1.0; A[0, 1] = alpha_b;
	A[1, 0] = alpha_1; A[1, 1] = 1.0; A[1, 2] = alpha_1;
	for i in range(2, nf-2):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;

	# Nagarajan et al boundary scheme
	A[nf-2, nf-3] = alpha_1; A[nf-2, nf-2] = 1.0; A[nf-2, nf-1] = alpha_1;
	A[nf-1, nf-2] = alpha_b; A[nf-1, nf-1] = 1.0;

	# build B
	# Nagarajan et al boundary scheme
	B[0,0] = a_b; B[0, 1] = b_b; B[0, 2] = c_b;
	B[1,0] = a_1/2.0; B[1, 1] = a_1/2.0;
	for i in range(2, nf-2):
		B[i, i-2] = b/2.0;
		B[i, i-1] = a/2.0;
		B[i, i] = a/2.0;
		B[i, i+1] = b/2.0;

	# Nagarajan et al boundary scheme
	B[nf-2, nc-2] = a_1/2.0; B[nf-2, nc-1] = a_1/2.0;
	B[nf-1,nc-3] = c_b; B[nf-1, nc-2] = b_b; B[nf-1, nc-1] = a_b;
	return A, B
#######################################

#######################################
# use a different function to create matrices
def create_matrices_interpolate_1_c2f(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 1st derivative
	from cell-center to cell-face (c2f).

    Inputs:
    x: grid (where f is defined), i.e., co-ordinates of cell-centers only

    Outputs:
    A, B to solve A varphi' = B psi
	"""
	nc = len(x)
	dx = (x[nc-1] - x[0])/(nc - 1)
	nf = nc + 1; # nface = nc + 1f

	# declare matrices
	A = np.zeros((nf, nf)); B = np.zeros((nf, nc));

    # create matrices
	alpha_b = 0.0; # alpha for boundary scheme Nagarajan et al.
	a_b = (-1.0/24.0)*(23.0*alpha_b + 71.0)
	b_b = (1.0/8.0)*(7.0*alpha_b + 47.0)
	c_b = (1.0/8.0)*(alpha_b - 31.0)
	d_b = (1.0/24.0)*(-alpha_b + 23.0)

	alpha_1 = 1.0/22.0    # alpha_1 for j = 1 and j = nc (second and second last to boundary)
	a_1 = (3.0/8.0)*(3.0 - 2.0*alpha_1)
	b_1 = 0.0

	alpha = 9.0/62.0    # alpha for inside the domain
	a = (3.0/8.0)*(3.0 - 2.0*alpha)
	b = (1.0/8.0)*(-1.0 + 22.0*alpha)

	# weights for discrete conservation
	w2 = 1.1 # stable free choice according to Nagarajan et. al.

	# solve the linear system given in Appendix B to find w0, w1, alpha_2
	# form the 3x3 system - M x = weights_rhs, where M = 3x3, weights_rhs is 3x1
	# and x = [w0, w1, alpha_2]'

	M = np.zeros((3, 3)); weights_rhs = np.zeros(3)

	M[0, 0] = b_b; M[0, 1] = a_1; M[0, 2] = (3.0/4.0)*w2;
	M[1, 0] = c_b; M[1, 1] = 0; M[1, 2] = (-3.0/4.0)*w2;
	M[2, 0] = d_b; M[2, 1] = 0; M[2, 2] = (11.0/12.0)*w2;

	weights_rhs[0] = (b/3.0) + (9.0*w2/8.0);
	weights_rhs[1] = a + (b/3.0) - (9.0*w2/8.0);
	weights_rhs[2] = (b/3.0)+(w2/24.0);

	temp_weights = solve(M, weights_rhs);
	w0 = temp_weights[0];
	w1 = temp_weights[1];
	alpha_2 = temp_weights[2];

	#print "w0 = ", w0
	#print "w1 = ", w1
	#print "alpha_2 = ", alpha_2

	a_2 = (3.0/8.0)*(3.0 - 2.0*alpha_2)
	b_2 = (1.0/8.0)*(-1.0 + 22.0*alpha_2)

	#print "a_2 = ", a_2
	#print "b_2 = ", b_2
	# build A
	# Nagarajan et al boundary scheme for discrete conservation
	A[0, 0] = 1.0*w0; A[0, 1] = alpha_b*w0;
	A[1, 0] = alpha_1*w1; A[1, 1] = 1.0*w1; A[1, 2] = alpha_1*w1;
	A[2, 1] = alpha_2*w2; A[2, 2] = 1.0*w2; A[2, 3] = alpha_2*w2;
	for i in range(3, nf-3):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;

	# Nagarajan et al boundary scheme for discrete conservation
	A[nf-3, nf-4] = alpha_2*w2; A[nf-3, nf-3] = 1.0*w2; A[nf-3, nf-2] = alpha_2*w2;
	A[nf-2, nf-3] = alpha_1*w1; A[nf-2, nf-2] = 1.0*w1; A[nf-2, nf-1] = alpha_1*w1;
	A[nf-1, nf-2] = alpha_b*w0; A[nf-1, nf-1] = 1.0*w0;

	# build B
	# Nagarajan et al boundary scheme
	B[0,0] = a_b*w0; B[0, 1] = b_b*w0; B[0, 2] = c_b*w0; B[0, 3] = d_b*w0;
	B[1,0] = -a_1*w1; B[1, 1] = a_1*w1;
	B[2,0] = w2*(-b_2/3.0); B[2, 1] = -w2*a_2; B[2, 2] = w2*a_2; B[2, 3] = w2*(b_2/3.0);

	for i in range(3, nf-3):
		B[i, i-2] = -b/3.0;
		B[i, i-1] = -a;
		B[i, i] = a;
		B[i, i+1] = b/3.0;

	# Nagarajan et al boundary scheme for discrete conservation
	B[nf-3, nc-4] = w2*(-b_2/3.0); B[nf-3, nc-3] = -w2*a_2; B[nf-3, nc-2] = w2*a_2; B[nf-3, nc-1] = w2*(b_2/3.0);
	B[nf-2, nc-2] = -a_1*w1; B[nf-2, nc-1] = a_1*w1;
	B[nf-1,nc-4] = -d_b*w0; B[nf-1,nc-3] = -c_b*w0; B[nf-1, nc-2] = -b_b*w0; B[nf-1, nc-1] = -a_b*w0;

	B = (1.0/dx)*B
	return A, B
#######################################
#######################################

def interpolate_c2f(f, x, A, B):

	"""
	Interpolate the 0th derivative (function value itself)
	from cell-center to cell-face.

	Inputs:
	f: discrete function
	x: grid (where f is defined), i.e., co-ordinates of cell-faces only
	A, B : obtained from create_matrices_interpolate_zero_f2c
	Outputs:
	f_at_face: interpolated f at cell-face
	"""

	nc = len(x)
	f_at_face	= np.zeros(nc+1)

	rhs = np.matmul(B, f)

	f_at_face = qrsolve(A, rhs)

	# f_at_center = solve(A, rhs)
	return f_at_face
#######################################
#######################################
# use a different function to create matrices
def create_matrices_interpolate_1_colloc(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 1st derivative
	on the collocated grid.

    Inputs:
    x: grid (where f is defined)

    Outputs:
    A, B to solve A varphi' = B psi
	"""
	n = len(x)
	dx = (x[n-1] - x[0])/(n - 1)

	# declare matrices
	A = np.zeros((n, n)); B = np.zeros((n, n));

    # create matrices
	alpha_b = 3.0; # alpha for boundary scheme Lele.
	a_b = (-1.0/6.0)*(11.0 + 2.0*alpha_b)
	b_b = (1.0/2.0)*(6.0 - alpha_b)
	c_b = (1.0/2.0)*(2.0*alpha_b - 3.0)
	d_b = (1.0/6.0)*(2.0 - alpha_b)

	# alpha_1 for j = 1 and j = n (second and second last to boundary)
	# Std 4th order Pade scheme
	alpha_1 = 1.0/4.0
	a_1 = (3.0/4.0)
	b_1 = 0.0

	alpha = 1.0/3.0    # alpha for inside the domain
	# qHat, rHat in Lele is a and b here
	a = (2.0/3.0)*(alpha + 2.0)*0.5
	b = (1.0/3.0)*(4.0*alpha - 1.0)*0.25

	# weights for discrete conservation
	# w0 here is w1 in Lele and wj here is w(j-1) in Lele
	# alpha_2 = alpha'' obtained in Lele
	alpha_2 =  ((40.0*alpha - 1.0)*b_b + 7.0*(4.0*alpha - 1.0)*d_b)\
	/ (16.0*(alpha + 2.0)*b_b + 8.0*(1.0 - 4.0*alpha)*d_b)

	w0 = (2.0*alpha + 1.0)/(2.0*(b_b + d_b))
	w1 = (8.0*alpha + 7.0)/9.0 - (6.0*(2.0*alpha + 1.0)*c_b)/(9.0*(b_b + d_b))
	w2 = (4.0*(alpha + 2.0)*b_b + 2.0*(1.0 - 4.0*alpha)*d_b)/(9.0*(b_b + d_b))

	#print "w0 = ", w0
	#print "w1 = ", w1
	#print "alpha_2 = ", alpha_2

	a_2 = (2.0/3.0)*(alpha_2 + 2.0)*0.5
	b_2 = (1.0/3.0)*(4.0*alpha_2 - 1.0)*0.25

	#print "a_2 = ", a_2
	#print "b_2 = ", b_2

	# build A
	# Lele boundary scheme for discrete conservation
	A[0, 0] = 1.0*w0; A[0, 1] = alpha_b*w0;
	A[1, 0] = alpha_1*w1; A[1, 1] = 1.0*w1; A[1, 2] = alpha_1*w1;
	A[2, 1] = alpha_2*w2; A[2, 2] = 1.0*w2; A[2, 3] = alpha_2*w2;
	for i in range(3, n-3):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;

	# Lele boundary scheme for discrete conservation
	A[n-3, n-4] = alpha_2*w2; A[n-3, n-3] = 1.0*w2; A[n-3, n-2] = alpha_2*w2;
	A[n-2, n-3] = alpha_1*w1; A[n-2, n-2] = 1.0*w1; A[n-2, n-1] = alpha_1*w1;
	A[n-1, n-2] = alpha_b*w0; A[n-1, n-1] = 1.0*w0;

	# build B
	# Lele boundary scheme
	B[0,0] = a_b*w0; B[0, 1] = b_b*w0; B[0, 2] = c_b*w0; B[0, 3] = d_b*w0;
	B[1,0] = -a_1*w1; B[1, 2] = a_1*w1;
	B[2,0] = -w2*b_2; B[2, 1] = -w2*a_2; B[2, 3] = w2*a_2; B[2, 4] = w2*(b_2);

	for i in range(3, n-3):
		B[i, i-2] = -b;
		B[i, i-1] = -a;
		B[i, i+1] = a;
		B[i, i+2] = b;

	# Lele boundary scheme for discrete conservation
	B[n-3, n-5] = -w2*(b_2); B[n-3, n-4] = -w2*a_2; B[n-3, n-2] = w2*a_2; B[n-3, n-1] = w2*(b_2);
	B[n-2, n-3] = -a_1*w1; B[n-2, n-1] = a_1*w1;
	B[n-1,n-4] = -d_b*w0; B[n-1,n-3] = -c_b*w0; B[n-1, n-2] = -b_b*w0; B[n-1, n-1] = -a_b*w0;

	B = (1.0/dx)*B
	return A, B
#######################################
#######################################
# use a different function to create matrices
def create_matrices_interpolate_2_colloc(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 2nd derivative
	on the collocated grid.

    Inputs:
    x: grid (where f is defined)

    Outputs:
    A, B to solve A varphi'' = B psi
	"""
	n = len(x)
	dx = (x[n-1] - x[0])/(n - 1)

	# declare matrices
	A = np.zeros((n, n)); B = np.zeros((n, n));

    # create matrices
	alpha_b = 10.0; # alpha for boundary scheme Lele.
	a_b = (1.0/12.0)*(11.0*alpha_b + 35.0)
	b_b = (-1.0/3.0)*(5.0*alpha_b + 26.0)
	c_b = (1.0/2.0)*(alpha_b + 19.0)
	d_b = (1.0/3.0)*(alpha_b - 14.0)
	e_b = (1.0/12.0)*(11.0 - alpha_b)

	# print "a_b = ", a_b
	# alpha_1 for j = 1 and j = n (second and second last to boundary)
	# Std 4th order Pade scheme
	alpha_1 = 1.0/10.0
	a_1 = (-12.0/5.0)
	b_1 = (6.0/5.0)

	alpha =2.0/11.0    # alpha for inside the domain
	# qHat, rHat in Lele is a and b here
	a = (4.0/3.0)*(1.0 - alpha)
	b = (1.0/3.0)*(-1 + 10.0*alpha)

	# weights for discrete conservation - to be updated
	w0 = 1.0;
	w1 = 1.0;
	w2 = 1.0;

	# build A
	# Lele boundary scheme for discrete conservation
	A[0, 0] = 1.0*w0; A[0, 1] = alpha_b*w0;
	A[1, 0] = alpha_1*w1; A[1, 1] = 1.0*w1; A[1, 2] = alpha_1*w1;
	A[2, 1] = alpha*w2; A[2, 2] = 1.0*w2; A[2, 3] = alpha*w2;
	for i in range(3, n-3):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;

	# Lele boundary scheme for discrete conservation
	A[n-3, n-4] = alpha*w2; A[n-3, n-3] = 1.0*w2; A[n-3, n-2] = alpha*w2;
	A[n-2, n-3] = alpha_1*w1; A[n-2, n-2] = 1.0*w1; A[n-2, n-1] = alpha_1*w1;
	A[n-1, n-2] = alpha_b*w0; A[n-1, n-1] = 1.0*w0;

	# build B
	# Lele boundary scheme
	B[0,0] = a_b*w0; B[0, 1] = b_b*w0; B[0, 2] = c_b*w0; B[0, 3] = d_b*w0; B[0, 4] = e_b*w0;
	B[1,0] = b_1*w1; B[1, 1] = a_1*w1; B[1, 2] = b_1*w1;
	B[2,0] = w2*b*0.25; B[2, 1] = w2*a; B[2, 2] = w2*(-2.0*a - 0.5*b); B[2, 3] = w2*a; B[2, 4] = w2*b*0.25;

	for i in range(3, n-3):
		B[i, i-2]	= 	0.25*b
		B[i, i-1]	= 	a
		B[i, i]		= 	-2.0*a - 0.5*b
		B[i, i+1]	= 	a
		B[i, i+2]	=	0.25*b

	# Lele boundary scheme for discrete conservation
	B[n-3, n-5] = w2*b*0.25; B[n-3, n-4] = w2*a; B[n-3, n-3] =  w2*(-2.0*a - 0.5*b) ;B[n-3, n-2] = w2*a; B[n-3, n-1] = w2*b*0.25;
	B[n-2, n-3] = b_1*w1; B[n-2, n-2] = a_1*w1; B[n-2, n-1] = b_1*w1;
	B[n-1, n-5] = e_b*w0; B[n-1, n-4] = d_b*w0; B[n-1, n-3] = c_b*w0; B[n-1, n-2] = b_b*w0; B[n-1, n-1] = a_b*w0;

	B = (1.0/(dx*dx))*B
	return A, B
#######################################
#######################################
# use a different function to create matrices
def create_matrices_interpolate_1_colloc_pbc(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 1st derivative
	on the collocated grid with periodic boundary conditions (pbc).

    Inputs:
    x: grid (where f is defined)

    Outputs:
    A, B to solve A varphi' = B varphi
	"""
	n = len(x)
	dx = (x[n-1] - x[0])/(n - 1)

	# declare matrices
	A = np.zeros((n, n)); B = np.zeros((n, n));

	# create matrices
	alpha = 1.0/3.0    # alpha for inside the domain
	# qHat, rHat in Lele is a and b here
	a = (2.0/3.0)*(alpha + 2.0)
	b = (1.0/3.0)*(4.0*alpha - 1.0)

	# build A
	# Lele boundary scheme for discrete conservation
	A[0, n-1] = alpha; A[0, 0] = 1.0; A[0, 1] = alpha;   # pbc
	for i in range(1, n-1):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;
	A[n-1, n-2] = alpha; A[n-1, n-1] = 1.0; A[n-1, 0] = alpha; # pbc

	# build B
	# pbc
	B[0, n-2] = -b*0.25; B[0, n-1] = -a*0.5; B[0, 1] = a*0.5; B[0, 2] = b*0.25;
	B[1, n-1] = -b*0.25; B[1, 0] = -a*0.5; B[1, 2] = a*0.5; B[1, 3] = b*0.25;

	for i in range(2, n-2):
		B[i, i-2] = -b*0.25;
		B[i, i-1] = -a*0.5;
		B[i, i+1] = a*0.5;
		B[i, i+2] = b*0.25;

	# pbc
	B[n-2, n-4] = -b*0.25; B[n-2, n-3] = -a*0.5; B[n-2, n-1] = a*0.5; B[n-2, 0] = b*0.25;
	B[n-1, n-3] = -b*0.25; B[n-1, n-2] = -a*0.5; B[n-1, 0] = a*0.5; B[n-1, 1] = b*0.25;

	B = (1.0/dx)*B
	return A, B
#######################################
#######################################
# use a different function to create matrices
def create_matrices_interpolate_2_colloc_pbc(x):
	"""
	The A and B matrices are constants and I need not define them
	again and again. Instead, I will just create them once and for all
	before the time loop and use them.

    Interpolate the 1st derivative
	on the collocated grid with periodic boundary conditions (pbc).

    Inputs:
    x: grid (where f is defined)

    Outputs:
    A, B to solve A varphi'' = B varphi
	"""
	n = len(x)
	dx = (x[n-1] - x[0])/(n - 1)

	# declare matrices
	A = np.zeros((n, n)); B = np.zeros((n, n));

    # create matrices
	alpha = 2.0/11.0    # alpha for inside the domain
	# qHat, rHat in Lele is a and b here
	a = (4.0/3.0)*(1.0 - alpha)
	b = (1.0/3.0)*(-1 + 10.0*alpha)

	# build A
	A[0, n-1] = alpha; A[0, 0] = 1.0; A[0, 1] = alpha;   # pbc
	for i in range(1, n-1):
		A[i, i-1]	= alpha;
		A[i, i]		= 1.0;
		A[i, i+1]	= alpha;
	A[n-1, n-2] = alpha; A[n-1, n-1] = 1.0; A[n-1, 0] = alpha; # pbc

	# build B
	# pbc
	B[0, n-2] = 0.25*b; B[0, n-1] = a; B[0, 0] = -2.0*a - 0.5*b; B[0, 1] = a; B[0, 2] = 0.25*b;
	B[1, n-1] = 0.25*b; B[1, 0] = a; B[1, 1] = -2.0*a - 0.5*b; B[1, 2] = a; B[1, 3] = 0.25*b;

	for i in range(2, n-2):
		B[i, i-2]	= 	0.25*b
		B[i, i-1]	= 	a
		B[i, i]		= 	-2.0*a - 0.5*b
		B[i, i+1]	= 	a
		B[i, i+2]	=	0.25*b

	# pbc
	B[n-2, n-4] = 0.25*b; B[n-2, n-3] = a; B[n-2, n-2] = -2.0*a - 0.5*b; B[n-2, n-1] = a; B[n-2, 0] = 0.25*b;
	B[n-1, n-3] = 0.25*b; B[n-1, n-2] = a; B[n-1, n-1] = -2.0*a - 0.5*b; B[n-1, 0] = a; B[n-1, 1] = 0.25*b;

	B = (1.0/(dx*dx))*B
	return A, B
#######################################

#######################################

def interpolate_colloc(f, x, A, B):

	"""
	Interpolate on the collocated grid

	Inputs:
	f: discrete function
	x: grid (where f is defined), i.e., co-ordinates of cell-faces only
	A, B : obtained from create_matrices_interpolate_zero_f2c
	Outputs:
	f_at_face: interpolated f at cell-face
	"""

	n = len(x)
	f_at_face	= np.zeros(n)

	rhs = np.matmul(B, f)

	# f_at_face = qrsolve(A, rhs)

	f_at_center = solve(A, rhs)
	return f_at_face
#######################################
