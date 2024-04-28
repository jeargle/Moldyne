Moldyne
=======

Molecular dynamics and monte carlo simulator written in julia.

This project is based off of lessons learned from *Understanding Molecular Dynamics* by Frenkel and Smit, *The Art of Molecular Dynamics Simulations* by Rapaport, and a Computational Statistical Mechanics course on coursera.


Input
------

The only input file is a PDB (Protein Data Bank) file with ATOM lines containing the atomic positions (2D or 3D).


Output
------

Two output text files are generated:

1. .out - summary stats per timestep
2. .xyz - coordinates per timestep; can be loaded into VMD as a trajectory


Dimension
---------

The simulation can run in 2D or 3D mode.  Both modes can use the same atomic position files, but 2D mode ignores the third position coordinate.


Dependencies
------------

* Distributions
* Plots
* Printf
* Random
* SpecialFunctions
