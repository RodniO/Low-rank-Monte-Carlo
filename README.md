# Projective volume low-rank approximation

Hello!
This repository contains Fortran code of the methods from
[2] A.I. Osinsky. Rectangular maximum volume and projective volume search algorithms // arXiv 1809.02334 (Submitted on 7 Sep 2018)

There are vector (ModVec) and matrix (ModMtrx) modules, which can be quite useful. However, they still contain a lot of trash (like old version of maxvol-related algorithms), so use them with caution.

The code of Dominant-R differs from the paper. I speed it up by not computing QR at all and using different recalculation formulas.

The use of all the new algorithms is illustrated in incfiles/Example.f90

To compile and run it, use either

./Launcher

(Requires gfortran, BLAS and LAPACK)

or

make

./Main.exe

(Requires ifort and mkl)

The output should look like this:

 Welcome to the low-rank approximation example!
 
 We are going to use different methods to
 compute fast low-rank approximations and compare them.
 If you get no error here, you are fine!
 
 We generate a 2000 by 2000 random matrix...
 Done!
 
 We set the first 10 singular values to be equal to 100 and others to 1.
 Done!
 
 We seek rank 10 approximation, so Frobenius norm error of SVD is
 SVD error:   44.6094160463909     
 
 Next we perform MAXVOL approximation
 MAXVOL time:  2.636848483234644E-003
 MAXVOL error:   69.8891690343534     
 
 Then let us try Householder-based MAXVOL2
 to construct FAST CGR.
 We add rows and columns up to 10*2 = 20
 Remember, that we search in rows and columns from MAXVOL
 FAST CGR time:  1.975158927962184E-003
 FAST CGR error:   54.7672233617187     
 
 Nobody needs maxvol-rect separately, so let us use MAXVOL-PROJ
 We discard previous rows and columns
 (To illustrate that maxvol-proj can work without initialization)
 And try to construct approximation from random start
 MAXVOL-PROJ time:  1.362081291154027E-002
 MAXVOL-PROJ error:   54.6947916025706     
 
 Finally, we construct Strong Rank Revealing QR with Dominant-R
 DOMINANT-R RRQR time:  0.130547259701416     
 DOMINANT-R RRQR error:   55.3540884942439 
