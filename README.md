# Projective volume low-rank approximation

Hello!
This repository contains a Fortran code of the low-rank Monte Carlo methods for (generalized) Smoluchowski equations from
[1] A.I. Osinsky. Low-rank Monte Carlo for Smoluchowski-class equations // Journal of Computational Physics 506(1), 112942 (2024).

and for Boltzmann equations from
[2] A.S. Bodrova, A.I. Osinsky. Anomalous diffusion in polydisperse granular gases: Monte Carlo simulations // arXiv 2403.13772 (Submitted on 20 Mar 2024).

Module ModMonte.mod contains all the corresponding subroutines. Other modules are from "projective volume low rank" repository (this repository is essentially a fork) and will be used, when ODE solvers for (generalized) Smoluchowski equations are added here.

The use of all of the simulation methods is illustrated in incfiles/ExampleM.f90.

To compile and run it, use either

make gnu_run

(Requires gfortran, BLAS and LAPACK)

or

make run

(Requires ifort and mkl)

Structure:

The code is organized to enable compilation of multiple programs with a single makefile. This is achieved through 'src/incfiles' directory, where one can write functions and then include and call them in 'Main.f90'. Makefile pays attention only to files with '.f90' extension, so remember to name them appropriately.
Module and object files are saved in 'obj_gnu' and 'obj_intel' folders (created during compilation), enabling to have gfortran and ifort versions compiled and saved simultaneously.

Large sizes:

Code uses 'allocatable' arrays, which are usually put in stack. Therefore, for large amounts of data, you may require to use 'ulimit -s unlimited' to allow unlimited stack size.

Debug:

debugrun.sh script uses gfortran with debug options to recompile and run everything. Remember to run 'make clean' after debug before running 'make gnu'.
