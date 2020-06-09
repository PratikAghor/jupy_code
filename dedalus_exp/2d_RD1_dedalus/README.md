# 2d reaction-diffusion equation in dedalus
## Author: Pratik Aghor


### Files and usage 
* ```RD_1.py``` uses periodic BC's in both x- and y-directions
* ```RD_2.py``` uses homogeneous Neumann BC's in both x- and y-directions (doesn't work as dedalus apparently cannot handle Chebyshev in two directions)
* ```RD_3.py``` uses periodic BC's in the x-direction and homogeneous Neumann BC's in the y-direction

* after activating dedalus type ```python3 RD_3.py```
* ```python3 plot_slices.py snapshots/snapshots_s20/snapshots_s20_p0.h5```
* Following seems to be a pretty robust steady state

![ss1]("write_001000.png")

* Need to work more on this. 
