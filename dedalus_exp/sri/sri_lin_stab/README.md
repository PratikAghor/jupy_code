# Linear Stability of Stratified Taylor-Couette Flow (TCF)
# Author: Pratik Aghor

* [Robins et al., 2020] Robins, J.M., Kersalé, E., and Jones, C. A. (2020). Viscous and Inviscid Strato-rotational Instability.
Journal of fluid mechanics.
* [Khorrami, 1990] A Cheb Spectral Collocation Method Using Staggered Grid for the Stability 
of Cylindrical Flows, Int Journal of Numerical Methods in Fluids.
* For details, read the blog. (do ```pdflatex blogname.tex``` in the terminal)

## Linear Stability of Stratified TCF:
* Started writing a modular python code 
* **NOTE**: for everything to work smoothly, change sys path from ```sys.path.append('/home/aghor/Aghor/UNH/independent/jupy_codes_local/dedalus_exp/sri/sri_lin_stab')``` to whatever the current path is. 
* To run tests, do ```bash build_tests.sh```
* To run the main code, do ```bash build.sh```

## cheb.py
* Given ```N```, returns the Gauss-Lobatto (GL) grid as well as the Chebyshev differentiation matrix
* Tested, works well. 

## interpolants.py
* The function ```gl2g(N)``` returns the matrix GL2G which interpolates a function at the GL-grid to the Gauss grid (G). 
* The function ```g2gl(N)``` returns the matrix G2GL which interpolates a function at the G-grid to the GL-grid. 
* We will need these interpolation matrices for the continuity equation. 
* As it has no boundary conditions, we will enforce the linearized conituity equation at the G-grid
* Then, using the G2GL matrix, we will obtain values at the GL-grid. 
* This way, we don't have to invent new boundary conditions (BC's) for density. 
* Tested both functions, work well for ```N > ~ 63```. For more oscillatory functions, more grid-points will be required. 
* Added a matrix ```D_g2gl(N)``` to obtain the first derivative at the GL-grid from a field defined at the G-grid. 
* **NOTE**: We prefer using square matrices, so the matrices are padded with zeros. Hence going from G2GL, we first must pad the vector with a zero at the end and then use matmul. Similarly, when going from GL2G, we must ignore the last element obtained from interpolations. 
