MB-pol is described here: <http://pubs.acs.org/doi/abs/10.1021/ct500079y>

To use MB-pol:

1. Go to the mbpol subdirectory and build `libmbpol.a` [adjust the `Makefile`
   as needed]; by default the `Makefile` uses intel compiler and compiles a
   thread safe version, this is going to take some time and memory and produce
   `libmbpol.a`

2. Fortran/C++ linkage examples are in the corresponding sub-folders
   (you may need to adjust the compiler flags); by default the `Makefile` uses intel 
   compiler.


There are minor differences between this version and the code published in 
the SI of <http://pubs.acs.org/doi/abs/10.1021/ct500079y>

1. If a monomer distortion energy is more than 100 kcal/mol, the 2B and higher
   many-body energy/gradient calculations are disabled, and the total energy
   becomes the sum of the monomer distortion energies. This is near the 
   homolytic dissociation energy of a (gas-phase) water molecule, and MB-pol
   is not meant to describe this region.

2. Minor performance bugfix in the calculation of the three-body terms. Bug did
   not affect results obtained with the previous code, however, the short-ranged
   three-body terms of excluded (long-range) trimers were needlessly calculated,
   slowing performance.
