
!Calculates matrix element for i-th row and j-th column.
pure function Exp_elem(i, j, param) Result(res)
  USE ModMtrx
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  
  res = exp(-0.01d0*i*j)
end

subroutine ExampleP()
  USE ModAppr
  Type(Mtrx) A, B, U, S, V, param, C, UR, E
  Type(IntVec) per1, per2
  Integer(4) n, k
  Double precision dsecnd, time
  
  !WELCOMING AND INITIALIZATION
  print *, ''
  print *, 'Welcome to the third example!'
  print *, ''
  print *, 'Now we construct nonnegative approximation of exp(-0.01*i*j).'
  print *, ''
  
  !Desired size
  n = 1000
  !Desired rank
  k = 10
  !We do not need parameters, but it is not good to leave it uninitialized
  call param%init(1,1)
  !We will also need permutations for Acols to work with
  call per1%perm(n)
  call per2%perm(n)
  
  print *, 'Singular values are exponentially decreasing'
  print *, 'so it is much better to use maxvol'
  print *, 'and then perform truncated svd of the approximation.'
  print *, ''
  
  !Just use all the columns to create the matrix.
  !We won't need it later, it is just used to see how SVD performs.
  A = Acols(Exp_elem, n, n, per1, per2, param)
  
  write(*,'(A,I0,A)') ' First of all let us look at the SVD error of rank ', k, ' approximation.'
  write(*,'(A,I0,A,I0,A)') ' Computing SVD of ', n ,' by ', n, ' matrix...'
  print *, ''
  
  call A%svd(U, S, V)
  !Only leave singular values after k
  S = S%subarray(n, n, k+1, k+1)
  print *, 'SVD error:', S%fnorm()
  B = U%subarray(n,k)*S%subarray(k,k)*V%subarray(k,n)
  print *, 'Minimum element:', minval(B%d)
  print *, ''
  
  print *, 'Next we perform alternating projections to make it nonnegative.'
  !Let's calculate the time
  time = dsecnd()
  
  !We use acc_type=0 to make it faster (it relies on premaxvol)
  call PositCUR(Exp_elem, param, n, n, k, 2*k, 3*k, 0.0d0, 0.0d0, 1.0d0, C, UR, 0, acc_type = 0)
  
  time = dsecnd() - time
  print *, 'Approximation time:', time
  
  B = C*UR
  E = A - B
  print *, 'Approximation error:', E%fnorm()
  print *, 'Minimum element is now', minval(B%d)
end
