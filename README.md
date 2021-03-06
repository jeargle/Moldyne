Moldyne
=======

Molecular dynamics program written in ruby.

It is based off of lessons learned in *Understanding Molecular Dynamics* by Frenkel and Smit as well as *The Art of Molecular Dynamics Simulations* by Rapaport.


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
