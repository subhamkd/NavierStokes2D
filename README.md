# NavierStokes2D
This repository contains a flow solver that solves the Navier Stokes equations for a 2D incompressible laminar flow in a rectangular domain. The solver was written using the "12-steps to Navier Stokes" tutorial by Professor Lorena Barba as given in Barba, Lorena A., and Forsyth, Gilbert F. (2018). CFD Python: the 12 steps to Navier-Stokes equations. Journal of Open Source Education, 1(9), 21, https://doi.org/10.21105/jose.00021

The functions that solve the pressure poisson equation and the momentum equations are defined in the solver.py script. 

The results are saved after a desired number of timesteps as .csv files in a new folder named '__Results__'.  
The saved results are used to make an animation of the transient behaviour, which is saved as a gif in the __Gifs__ folder.
