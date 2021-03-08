import numpy as np
'''
Inviscid SRI linear stability
author: Pratik Aghor
Ref: Viscous and inviscid strato-rotational instabilities, Robins et. al., JFM 2020.
'''
###########################################
def spurious_function(xs):
    '''
    Determine if the given eigenfunction is spurious or not.
    big_sum = the sum of squares of absolute values of the coeffs in the first half of the Cheb expansion.
    small_sum = the sum of squares of absolute values of the coeffs in the second half of the Cheb expansion.

    if small_sum/big_sum > limit, then the eigenfunction is spurious. Filter it out.
    '''
    n = len(xs)
    half_len = int(n/2.)

    big_sum = 0.0 # first half
    small_sum = 0.0 # second half
    for j in range(0:half_len):
        big_sum = big_sum + (abs(xs[j]))**2
        small_sum = small_sum + (abs(xs[half_len + j]))**2
    ratio = small_sum/big_sum

    return ratio
###########################################
def spurious_filter(eigvals_, eigvecs_):
    """
    filter sputious eigenstuff
    safety = 1 => safe
    safety = 0 => spurious
    """
    limit = 1e-4
    nr = int(len(eigvals_)/5) # nr = Nr +1
    safety = np.ones(len(eigvals_)); # start with all safe
    for j in range(0, len(eigvals_)):
        u_ = eigvecs_[0:nr, j]
        v_ = eigvecs_[nr:2*nr, j]
        w_ = eigvecs_[2*nr:3*nr, j]
        rho_ = eigvecs_[3*nr:4*nr, j]
        p_ = eigvecs_[4*nr:5*nr, j]

        u_ratio = spurious_function(u_)
        v_ratio = spurious_function(v_)
        w_ratio = spurious_function(w_)
        rho_ratio = spurious_function(rho_)
        p_ratio = spurious_function(p_)

        # if even one condition fails, flag spurious
        if(u_ratio > limit or v_ratio > limit or w_ratio > limit or rho_ratio > limit or p_ratio > limit):
            safety[j] = 0;

    """
    now filter out spurious eigenstuff
    """
    return safety
###########################################
