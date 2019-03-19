
subroutine Example()
  USE ModMtrx
  Type(Mtrx) U, V, S, A, C, CA, AR, ART, E, Ahat, A1, A2, AB, RT, R, Q
  Type(Vector) pert1, pert2, pert3, perti
  Integer(4) n, k, maxswaps, swapsmade, i, j
  Double precision f, dsecnd, time
  
  !WELCOMING AND INITIALIZATION
  print *, 'Welcome to the low-rank approximation example!'
  print *, ''
  print *, 'We are going to use different methods to'
  print *, 'compute fast low-rank approximations and compare them.'
  print *, 'If you get no error here, you are fine!'
  print *, ''
  print *, 'We generate a 1000 by 1000 random matrix...'
  n = 1000
  call U%random(n)
  call V%random(n)
  call S%init(n,n)
  print *, 'Done!'
  print *, ''
  print *, 'We set the first 10 singular values to be equal to 100 and others to 1.'
  k = 10
  do i = 1, k
    S%d(i,i) = 100.0d0
  end do
  do i = k+1, n
    S%d(i,i) = 1.0d0
  end do
  A = U*S*V
  print *, 'Done!'
  print *, ''
  print *, 'We seek rank 10 approximation, so Frobenius norm error of SVD is'
  f = sqrt(dble(n-k))
  print *, 'SVD error:', f
  print *, ''
  
  !MAXVOL
  print *, 'Next we perform MAXVOL approximation'
  !We initialize the pert1 and pert2, which store permutations of rows and columns
  !We apply this permutations to put the maximum volume submatrix in the top left corner.
  !Some algorithms do it automatically.
  !You can uncomment some alternatives and comment the other: they do the same.
  call pert1%pert(n)
  call pert2%pert(n)
  !We will also need identity permutation
  call perti%pert(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  !Here we use 'cmaxvol', which is maxvol in columns
  maxswaps = 2*k
  !Let's do 4 steps
  do i = 1, 2
    !Select first k columns
    !Instead of Arows and Acols one can use any other
    !functions, which read or calculate values of A only
    !when called and do not use any precalculated information.
    !So A is not needed to be stored.
    C = Acols(A, k, pert1, pert2)
    !Use column version of maxvol.
    !One can also get the number of swaps made: swapsmade
    !One can use swapsmade=0 as a stopping criterion.
    !And C \hat A^{-1}: CA
    !Hereinafter all parameters after sizes are optional
    call C%cmaxvol(pert1, swapsmade, maxswaps, CA)
    !Select and transpose first k rows
    RT = .T.Arows(A, k, pert1, pert2)
    !We again use column version and swap columns
    !That's why we have 2 instead of 1
    call RT%cmaxvol(pert2, swapsmade, maxswaps, ART)
  end do
  !Rows should coincide with the rows of A
  C = Acols(A, k, perti, pert2)
  AR = .T.ART
  !We swap the columns back to make them coincide with the columns of A
  call AR%pertcols(pert2, 2)
  !Now our low-rank approximation is C \hat A^{-1} R = C*AR
  
  time = dsecnd() - time
  print *, 'MAXVOL time:', time
  
  E = A - C*AR
  print *, 'MAXVOL error:', E%fnorm()
  print *, ''
  
  time = dsecnd()
  
  !MAXVOL2
  print *, 'Then let us try Householder-based MAXVOL2'
  print *, 'to construct FAST CGR.'
  print *, 'We add rows and columns up to 10*2 = 20'
  print *, 'Remember, that we search in rows and columns from MAXVOL'
  !We need to start from good rows and columns, so
  !that the top left k by k matrix has maximum volume.
  !Here we use already calculated C and A^{-1} R
  !We return back the permutations
  call C%pertrows(pert1, 1)
  call AR%pertcols(pert2, 1)
  !Can be applied to the entire matrix like A%hmaxvol2 too, if necessary.
  call C%hmaxvol2(1, k, 2*k, pert1, 0)
  !2 indicates that we swap columns (1 if rows)
  !1 indicates that our rows are already multiplied be A^-1 (0 if not)
  call AR%hmaxvol2(2, k, 2*k, pert2, 1)
  !We now need to use more rows and columns
  C = Acols(A, 2*k, perti, pert2)
  R = Arows(A, 2*k, pert1, pert2)
  !Moreover, we need to use PROJECTIVE VOLUME
  Ahat = R%subarray(2*k,2*k)
  call Ahat%svd(U, S, V)
  !We use k-pseudoinverse
  U = U%subarray(2*k,k)
  S = S%subarray(k,k)
  V = V%subarray(k,2*k)
  do i = 1, k
    S%d(i,i) = 1.0d0/S%d(i,i)
  end do
  call R%pertcols(pert2, 2)
  !Multiplication in stable order
  AR = ((.T.V)*S) * ((.T.U)*R)
  
  time = dsecnd() - time
  print *, 'FAST CGR time:', time
  
  E = A - C*AR
  print *, 'FAST CGR error:', E%fnorm()
  print *, ''
  
  !MAXVOL-PROJ
  print *, 'Nobody needs maxvol-rect separately, so let us use MAXVOL-PROJ'
  print *, 'We discard previous rows and columns'
  print *, 'And try to construct approximation from the beginning'
  call pert1%deinit()
  call pert2%deinit()
  call pert1%pert(n)
  call pert2%pert(n)
  !We also need to keep some trash
  call pert3%pert(n)
  
  time = dsecnd()
  
  !We do essentially the same we have been doing with 'cmaxvol'
  !so no comment

  !maxvol-rect of size 2k x k
  do i = 1, 2
    R = Arows(A, 2*k,pert1, pert3)
    call R%dominantr(2, k, 2*k, pert3, swapsmade)
    C = Acols(A, k, pert1, pert3)
    !1 for swaps of rows; 0 for 0 rows kept unswapped
    call C%dominantc(1, k, 2*k, 0, pert1, pert3, swapsmade)
  end do
  call pert3%deinit()
  call pert3%pert(n)
  !maxvol-rect of size k x 2k
  do i = 1, 2
    C = Acols(A, 2*k, pert3, pert2)
    call C%dominantr(1, k, 2*k, pert3, swapsmade)
    R = Arows(A, k, pert3, pert2)
    call R%dominantc(2, k, 2*k, 0, pert3, pert2, swapsmade)
  end do
  call pert3%deinit()
  !Projective volume business again
  C = Acols(A, 2*k, perti, pert2)
  R = Arows(A, 2*k, pert1, pert2)
  Ahat = R%subarray(2*k,2*k)
  call Ahat%svd(U, S, V)
  U = U%subarray(2*k,k)
  S = S%subarray(k,k)
  V = V%subarray(k,2*k)
  do i = 1, k
    S%d(i,i) = 1.0d0/S%d(i,i)
  end do
  call R%pertcols(pert2, 2)
  AR = ((.T.V)*S) * ((.T.U)*R)
  
  time = dsecnd() - time
  print *, 'MAXVOL-PROJ time:', time
  
  E = A - C*AR
  print *, 'MAXVOL-PROJ error:', E%fnorm()
  print *, ''
  
  print *, 'Finally, we construct Strong Rank Revealing QR with Dominant-R'
  !Let's reinitialize column permutations
  call pert2%deinit()
  call pert2%pert(n)
  call AR%deinit()
  call AR%init(k, n)
  
  time = dsecnd()
  
  !We will work in the copy of A
  call A1%copy(A)
  !Let's limit ourselves to 2*k swaps
  call A1%dominantr(2, k, n, pert2, swapsmade, 2*k, Ahat, AB)
  
  Ahat = .T.Ahat
  !We do not need Q; pert3 and Q are dummy variables here
  call Ahat%halfLQ(R,pert3,Q)
  R = .T.R
  C = Acols(A, k, perti, pert2)
  Q = C*R
  AR%d(:,1:k) = 0
  do i = 1, k
    AR%d(i,i) = 1.0d0
  end do
  AR%d(:,k+1:n) = AB%d(:,1:n-k)
  AR = R%rtsolve(AR)
  call AR%pertcols(pert2, 2)
  
  time = dsecnd() - time
  print *, 'DOMINANT-R RRQR time:', time
  
  E = A - Q*AR
  print *, 'DOMINANT-R RRQR error:', E%fnorm()
end

function Arows(A, r, pert1, pert2) Result(res)
  USE ModMtrx
  Type(Mtrx), intent(in) :: A
  Integer(4), intent(in) :: r
  Type(Vector), intent(in) :: pert1, pert2
  Type(Mtrx) :: res
  Integer(4) i, j
  call res%init(r, A%m)
  do i = 1, r
    do j = 1, A%m
      res%d(i,j) = A%d(floor(pert1%d(i)+0.5d0), floor(pert2%d(j)+0.5d0))
    end do
  end do
end

function Acols(A, r, pert1, pert2) Result(res)
  USE ModMtrx
  Type(Mtrx), intent(in) :: A
  Integer(4), intent(in) :: r
  Type(Vector), intent(in) :: pert1, pert2
  Type(Mtrx) :: res
  Integer(4) i, j
  call res%init(A%n, r)
  do j = 1, r
    do i = 1, A%n
      res%d(i,j) = A%d(floor(pert1%d(i)+0.5d0), floor(pert2%d(j)+0.5d0))
    end do
  end do
end
