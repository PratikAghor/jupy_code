## Strato-Rotational Instability (SRI) in Dedalus:

### Pratik Aghor

* The plan is as follows: implement the linear stability analysis from `The Boussinesq approximation in rapidly rotating flows' by Lopez et. al. using dedalus's EVP solver
* Wrote a basic version, does not work yet. 

## sri_test_001:
* Go to dedalus mode.
* ```python3 sri_eigenvalue_3d_take_2.py``` gives the following result:
```
2020-07-28 19:44:16,990 __main__ 0/1 INFO :: Re:5.182e+03, eta:6.5000e-01, G:5.100e+01, epsilon:0.000e+00, Pr:4.3500e+00, mu:0.0000e+00
2020-07-28 19:44:16,990 __main__ 0/1 INFO :: Lz set to 2.007407e+00
2020-07-28 19:44:16,992 problems 0/1 INFO :: Solving EVP with homogeneity tolerance of 1.000e-10
2020-07-28 19:44:17,711 problems 0/1 INFO :: Solving EVP with homogeneity tolerance of 1.000e-10
/home/aghor/soft/dedalus_helpers/eigentools/eigentools/eigenproblem.py:232: ComplexWarning: Casting complex values to real discards the imaginary part
  indx = lambda1_and_indx[:, 1].astype(np.int)
2020-07-28 19:44:18,962 __main__ 0/1 INFO :: Growth rate = 1.347734575210464e+03; frequency = -1.314330107033370e+04

```

* TODO: verify the code with the published literature
* TODO: add functionality so that parameters can be passed as args

## sri_eigenvalue_3d_take_3:
* is a preliminary version of the SRI linear stability problem derived in `Viscous and inviscid strato-rotational instability by Robins et. al. (2020, JFM)'
* This does not consider the effects of Boussinesq approximation due to centrifugal forces and is in a regime where density variation is only considered in the axial direction. 
* Need to verify the results of Robins et. al. (2020, JFM) using a better version of this code.
