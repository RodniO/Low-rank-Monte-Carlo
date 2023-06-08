
!Calculates matrix element for i-th row and j-th column.
pure function Aelem(i, j, param) Result(res)
  USE ModMtrx
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  !Just uses parameters array as the input matrix
  res = param%d(i,j)
end

subroutine ExampleA()
  USE ModAppr
  Type(Mtrx) U, V, A, C, CA, AR, E, Ahat, A1, AB, ABin, R, Q
  Type(IntVec) per1, per2, per3, peri
  Type(Vector) cin, S
  Integer(4) n, k, maxsteps, maxswaps, swapsmade, i
  Double precision dsecnd, time
  
  !WELCOMING AND INITIALIZATION
  print *, ''
  print *, 'Welcome to the low-rank approximation example!'
  print *, ''
  print *, 'We are going to use different methods to'
  print *, 'compute fast low-rank approximations and compare them.'
  print *, 'If you get no error here, you are fine!'
  print *, ''
  
  !Desired size
  n = 1000
  !Desired rank
  k = 10
  !Maximum number of steps for maxvol and maxvol-rect
  maxsteps = 2
  !Maximum number of row and column swaps for maxvol and maxvol-rect
  maxswaps = 4*k
  
  write(*,'(A,I0,A,I0,A)') ' We create ', n ,' by ', n, ' matrix'
  print *, 'with random singular vectors (now generating)...'
  call U%random(n)
  call V%random(n)
  call S%init(n)
  print *, 'Done!'
  print *, ''
  
  write(*,'(A,I0,A)') ' We set the first ', k, ' singular values'
  print *, 'to 100 and others to 1...'
  do i = 1, k
    S%d(i) = 100.0d0
  end do
  do i = k+1, n
    S%d(i) = 1.0d0
  end do
  A = U*(S .dot. V)
  print *, 'Done!'
  print *, ''
  
  write(*,'(A,I0,A)') ' We seek rank ', k, ' approximation, so Frobenius norm error of SVD is'
  print *, 'SVD error:', sqrt(dble(n-k))
  print *, ''
  
  !MAXVOL
  print *, 'Next we perform MAXVOL approximation'
  !We initialize the perm1 and perm2, which store permutations of rows and columns
  !We apply this permutations to put the maximum volume submatrix in the top left corner.
  !Some algorithms do it automatically.
  call per1%perm(n)
  call per2%perm(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  !Find dominant k by k submatrix \hat A and construct C \hat A^{-1} R approximation
  call maxvol(Aelem, n, n, k, per1, per2, A, C, AR, maxsteps, maxswaps, CA, .true.)
  
  time = dsecnd() - time
  print *, 'MAXVOL time:', time
  
  E = A - C*AR
  print *, 'MAXVOL error:', E%fnorm()
  print *, ''
  
  time = dsecnd()
  
  !MAXVOL2
  print *, 'Then let us try Householder-based MAXVOL2'
  print *, 'to construct fast CUR.'
  write(*,'(A,I0,A,I0)') ' We add rows and columns up to ', k, '*2 = ', k*2
  print *, 'Remember, that we search in rows and columns from MAXVOL'
  
  !Increase submatrix size to 2k by 2k
  call maxvol2(Aelem, 2*k, 2*k, per1, per2, A, CA, AR, .true.)
  
  time = dsecnd() - time
  print *, 'FAST CUR time:', time
  
  E = A - CA*AR
  print *, 'FAST CUR error:', E%fnorm()
  print *, ''
  
  !MAXVOL-PROJ
  print *, 'Nobody needs maxvol-rect separately, so let us use MAXVOL-PROJ'
  print *, 'We discard previous rows and columns'
  print *, '(To illustrate that maxvol-proj can work without initialization)'
  print *, 'And try to construct approximation from random start'
  call per1%perm(n)
  call per2%perm(n)
  
  time = dsecnd()
  
  !Large projective volume submatrix search
  call maxvolproj(Aelem, n, n, k, 2*k, 2*k, per1, per2, A, C, AR, maxsteps, maxswaps)
  
  time = dsecnd() - time
  print *, 'MAXVOL-PROJ time:', time
  
  E = A - C*AR
  print *, 'MAXVOL-PROJ error:', E%fnorm()
  print *, ''
  
  !RRQR with DOMINANT_R
  print *, 'Finally, we construct Strong Rank Revealing QR with Dominant-R'
  !Let's reinitialize column permutations
  call per2%perm(n)
  !and rows
  call AR%init(k, n)
  
  time = dsecnd()
  
  !We'll need identity permuatation
  call peri%perm(n)
  !And some dummy permutation
  call per3%perm(n)
  
  !We will work in the copy of A
  call A1%copy(A)
  !We will use pre-maxvol to decrease time and increase accuracy
  !ABin and cin will be passed to Dominant-R to remove initialization of AB and c
  call A1%premaxvol(k, per2, ABin, cin)
  !Let's limit ourselves to 2*k swaps
  call A1%dominantr(2, k, n, per2, swapsmade, 2*k, Ahat, AB, ABin, cin)
  
  Ahat = .T.Ahat
  !We do not need Q; cin and Q are dummy variables here
  call Ahat%halfLQ(R, cin, Q)
  R = .T.R
  C = Acols(Aelem, N, k, peri, per2, A)
  Q = C*R
  do i = 1, k
    AR%d(i,i) = 1.0d0
  end do
  AR%d(:,k+1:n) = AB%d(:,1:n-k)
  AR = R%rtsolve(AR)
  call AR%permcols(per2, 2)
  
  time = dsecnd() - time
  print *, 'DOMINANT-R RRQR time:', time
  
  E = A - Q*AR
  print *, 'DOMINANT-R RRQR error:', E%fnorm()
end
