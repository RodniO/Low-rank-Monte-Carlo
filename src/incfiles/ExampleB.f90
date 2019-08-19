
!Calculates matrix element for i-th row and j-th column.
function Belem(i, j, param) Result(res)
  USE ModMtrx
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  !This is ballistic coagulation kernel.
  res = (i**(1.0d0/3.0d0) + j**(1.0d0/3.0d0))**2 * sqrt(1.0d0/i + 1.0d0/j)
end

subroutine ExampleB()
  USE ModAppr
  Type(Mtrx) A, U, S, V, param, C, UR, E
  Type(IntVec) per1, per2
  Integer(4) n, k, maxsteps, maxswaps
  Double precision dsecnd, time
  
  !WELCOMING AND INITIALIZATION
  print *, ''
  print *, 'Welcome to the second example!'
  print *, ''
  print *, 'Now we are to construct approximation of the ballistic coagulation kernel.'
  print *, ''
  
  !Desired size
  n = 1000
  !Desired rank
  k = 10
  !Maximum number of steps for maxvol
  maxsteps = 2
  !Maximum number of row and column swaps for maxvol
  maxswaps = 4*k
  !We do not need parameters, but it is not good to leave it uninitialized
  call param%init(1,1)
  !We will also need permutations for maxvol to work with
  call per1%perm(n)
  call per2%perm(n)
  
  print *, 'Singular values are exponentially decreasing'
  print *, 'so it is much better to use maxvol'
  print *, 'and then perform truncated svd of the approximation.'
  print *, ''
  
  !Just use all the columns to create the matrix.
  !We won't need it later, it is just used to see how SVD performs.
  A = Acols(Belem, n, n, per1, per2, param)
  
  print *, 'First of all let us look at the SVD error.'
  write(*,'(A,I0,A,I0,A)') ' Computing SVD of ', n ,' by ', n, ' matrix...'
  print *, ''
  
  call A%svd(U, S, V)
  !Only leave singular values after k
  S = S%subarray(n, n, k+1, k+1)
  print *, 'SVD error:', S%fnorm()
  print *, ''
  
  write(*,'(A,I0)') ' Next we perform MAXVOL approximation of rank ', k+5
  !Let's calculate the time
  time = dsecnd()
  
  !MAXVOL
  !Find dominant submatrix \hat A and construct C \hat A^{-1} R approximation
  call maxvol(Belem, n, n, k+5, per1, per2, param, C, UR, maxsteps, maxswaps)
  
  write(*,'(A,I0)') ' and use truncated SVD to decrease rank to ', k
  
  !TRUNCATED SVD
  call TruncateCUR(C, UR, k)
  
  time = dsecnd() - time
  print *, 'MAXVOL with TSVD time:', time
  
  E = A - C*UR
  print *, 'MAXVOL with TSVD error:', E%fnorm()
end
