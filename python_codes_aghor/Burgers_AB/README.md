# 1D Heat Equation Crank Nicholson
## Author: Pratik Aghor

* Time marching 1d Burger's equation using AB2 and AB4
* In order to run the codes, please do ```bash build.sh```

### Files and usage 

* ```params.py``` defines the parameters such as the grid and time step, etc.
* ```Burger_AB_functions.py``` defines the important matrices as well as boundary conditions. It also contains a code to solve Ax = b using QR factorization.
* ```Burger_AB2_main.py``` is the main file that does the AB2 time marching and saves the data in the data folder.
* ```Burger_AB4_main.py``` is the main file that does the AB4 time marching and saves the data in the data folder.
* ```post_process.py``` does post processing - reading and plotting the data, saving figures, etc.

