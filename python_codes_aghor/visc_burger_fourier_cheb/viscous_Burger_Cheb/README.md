# 1d Viscous Burger's Equation using Chebyshev Spectral Method
## Author: Pratik Aghor

* Time marching 1d viscous Burger's equation using Chebyshev spectral method
* In order to run the codes, please do ```bash build.sh```
* In order to run the test for ```cheb.py```, please do ```bash build_test.sh```

### Files and usage 

* ```params.py``` defines the parameters such as the grid and time step, etc.
* ```cheb.py``` contains the code to obtain the d()/dx operator (D-matrix) on the computational grid (Chebyshev-Gauss-Lobatto) and the Chebyshev-Gauss-Lobatto grid.
* ```test_cheb.py``` has a test to validate ```cheb.py```
* Exact solution vs numerical solution
![cheb_test](test_cheb.png)
* ```qrsolve.py``` contains the code to solve Ax = b using QR factorization.
* ```burger_cheb_PratikAghor.py``` is the main file that does the time marching and saves the data in the data folder.
* ```post_process.py``` does post processing - reading and plotting the data.
![ut](ut_cheb.png)

