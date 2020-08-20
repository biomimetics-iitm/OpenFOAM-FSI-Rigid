# About:

myRBF is a massively parallel implementation of RBFMotionSolver. Here, one need not manually choose to decompose the domain to include the moving boundary control points in each decomposed domain region. 

# Instructions to follow:
Just merge the src folder here into the foam-extend-4.0/src folder. 
Now, run the following to compile the mesh motion solver:
1. fe40 (or the corresponding environment initialisation alias for Foam extend 4.0)
2. cd /src/dynamicMesh
3. AllClean
4. Allwmake

Voila! You are done including the modified RBF mesh motion solver. Now, Enjoy running FSI/Rigid Body Motion cases with this solver.  

