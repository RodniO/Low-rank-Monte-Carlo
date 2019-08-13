# Projective volume low-rank approximation

Hello!
This repository contains Fortran code of the methods from
[2] A.I. Osinsky. Rectangular maximum volume and projective volume search algorithms // arXiv 1809.02334 (Submitted on 7 Sep 2018)

Module ModAppr contains subroutines for low-rank approximations. There are also vector (ModVec) and matrix (ModMtrx) modules, which can be quite useful.

The use of all of the approximation algorithms is illustrated in incfiles/ExampleA.f90 and incfiles/ExampleB.f90

To compile and run it, use either

make gnu_run

(Requires gfortran, BLAS and LAPACK)

or

make run

(Requires ifort and mkl)

Structure:

The code is organised to enable compilatation of multiple programs with a single makefile. This is achieved through 'src/incfiles' directory, where one can write functions and then include and call them in 'Main.f90'. Makefile pays attention only to files with '.f90' extension, so remember to name them properly.
Module and object files are saved in 'obj_gnu' and 'obj_intel' folders (created during compilation), enabling to have gfortran and ifort versions compiled and saved simulaneously.

Large sizes:

Code uses 'allocatable' arrays, which are usually put in stack. Therefore, for large amounts of data you may require to use 'ulimit -s unlimited' to allow unlimited stack size.

Debug:

debugrun.sh script uses gfortran with debug options to recompile and run everything. Remember to run 'make clean' after debug before running 'make gnu'.
