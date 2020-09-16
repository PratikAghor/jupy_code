from numpy import pi,cos,arange, ones
import numpy as np

def g2gl(N):
    '''Interpolate field from Gauss to Gauss-Lobatto points
       Ref.: A Cheb Spectral Collocation Method Using Staggered Grid for
	   the Stability of Cylindrical Flows, Khorrami 1990.
       author: Pratik Aghor
    '''
    xc  = cos(pi*arange(0,N+1)/N) # xc = Chebyshev-Gauss-Lobatto grid
    xg  = cos(pi*arange(1, (2*N), 2)/(2.0*(N)))  # xg = Chebyshev-Gauss grid

    G2GL   = np.zeros((N+1, N+1))

    for j in range(0, N+1):
        for k in range(0, N):
            G2GL[j, k] = ((-1)**(j + k) * np.sqrt(1.0 - xg[k]*xg[k]))/(N*(xc[j] - xg[k]))

    return G2GL

def gl2g(N):
    '''Interpolate field from Gauss-Lobatto to Gauss points
       Ref.: A Cheb Spectral Collocation Method Using Staggered Grid for
	   the Stability of Cylindrical Flows, Khorrami 1990.
       author: Pratik Aghor
    '''
    xc  = cos(pi*arange(0,N+1)/N) # xc = Chebyshev-Gauss-Lobatto grid
    xg  = cos(pi*arange(1, (2*N), 2)/(2.0*(N)))  # xg = Chebyshev-Gauss grid
    c   = ones(N+1); c[0] = 2.0; c[N] = 2.0
    GL2G   = np.zeros((N+1, N+1))

    for j in range(0, N):
        for k in range(0, N+1):
            GL2G[j, k] = ((-1)**(j + k + 1) * np.sqrt(1.0 - xg[j]*xg[j]))/(c[k]*N*(xg[j] - xc[k]))

    return GL2G

def D_g2gl(N):
    '''Ontain first derivative at Gauss-Lobatto points from a field at Gauss points
       Ref.: A Cheb Spectral Collocation Method Using Staggered Grid for
	   the Stability of Cylindrical Flows, Khorrami 1990.
       author: Pratik Aghor
    '''
    xc  = cos(pi*arange(0,N+1)/N) # xc = Chebyshev-Gauss-Lobatto grid
    xg  = cos(pi*arange(1, (2*N), 2)/(2.0*(N)))  # xg = Chebyshev-Gauss grid

    E   = np.zeros((N+1, N+1))
    for k in range(0, N):
            E[0, k] = (((-1)**(k) * np.sqrt(1.0 - xg[k]*xg[k]))/(1.0*N))* \
            ( ( 1.0*N*N/(1.0 - xg[k]) ) - ( 1.0/(1.0 - xg[k])**2) )
            E[N, k] = (((-1)**(k+N) * np.sqrt(1.0 - xg[k]*xg[k]))/(1.0*N))* \
            ( ( 1.0*N*N/(1.0 + xg[k]) ) - ( 1.0/(1.0 + xg[k])**2) )

    for j in range(1, N):
        for k in range(0, N):
            E[j, k] = ((-1)**(j + k + 1) * np.sqrt(1.0 - xg[k]*xg[k]))/(N*(xc[j] - xg[k])**2)

    return E
