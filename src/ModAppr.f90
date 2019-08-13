Module ModAppr
  USE ModMtrx
  
  !Module for constructing fast approximations
  !Contains maxvol, maxvol2, maxvolproj and TruncateCUR

  abstract interface
    !Interface for function, returning matrix elements
    function elem_fun(i, j, param) Result(res)
      USE ModMtrx
      Integer(4), intent(in) :: i, j
      Type(Mtrx), intent(in) :: param
      Double precision :: res
    end function
  end interface
  
  Contains
  
!Returns r columns of length N with elements from Afun
function Acols(Afun, N, r, per1, per2, param) Result(res)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: N, r !Length and number of columns
  Type(IntVec), intent(in) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx) :: res !Output matrix of columns
  
  Integer(4) i, j
  
  call res%init(N, r)
  do j = 1, r
    do i = 1, N
      res%d(i,j) = Afun(per1%d(i), per2%d(j), param)
    end do
  end do
end
  
!Returns r rows of length N with elements from Afun
function Arows(Afun, r, N, per1, per2, param) Result(res)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: r, N !Number of rows and their length
  Type(IntVec), intent(in) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx) :: res !Output matrix of rows
  
  Integer(4) i, j
  
  call res%init(r, N)
  do i = 1, r
    do j = 1, N
      res%d(i,j) = Afun(per1%d(i), per2%d(j), param)
    end do
  end do
end

!Returns transposed matrix of r rows of length N with elements from Afun
function Arowst(Afun, r, N, per1, per2, param) Result(res)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: r, N !Number of rows and their length
  Type(IntVec), intent(in) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx) :: res !Output matrix of rows
  
  Integer(4) i, j
  
  call res%init(N, r)
  do j = 1, r
    do i = 1, N
      res%d(i,j) = Afun(per1%d(j), per2%d(i), param)
    end do
  end do
end
  
!Simplest CUR approximation with U = \hat A^{-1}
subroutine maxvol(Afun, M, N, rank, per1, per2, param, C, UR, maxsteps, maxswaps)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: M, N, rank !Sizes of A and the desired rank
  Type(IntVec) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx), intent(out) :: C !Columns for CUR approximation
  Type(Mtrx), intent(out) :: UR !UR of CUR approximation
  Integer(4), intent(in) :: maxsteps, maxswaps !Maximum number of steps and swaps

  Type(IntVec) peri !Identity permutation
  Type(Mtrx) RT !Transposed rows
  Type(Mtrx) URT !Transposed UR
  Integer(4) swapsmade1, swapsmade2 !Number of swaps in rows and columns
  Integer(4) i
  
  !Initialize identity permutation
  call peri%perm(M)
  
  !Here we use 'cmaxvol', which is maxvol in columns
  do i = 1, maxsteps
    !Select first r columns
    !Instead of Aelem one can use any other
    !function, which reads or calculates values of A only
    !when called and does not use any precalculated information.
    !So A is not needed to be stored.
    !Same goes if one wants to modify Arows or Acols.
    C = Acols(Afun, M, rank, per1, per2, param)
    !We use column version of maxvol.
    call C%cmaxvol(per1, swapsmade1, maxswaps)
    
    !Select first r rows (transposed)
    RT = Arowst(Afun, rank, N, per1, per2, param)
    !We again use column version and swap columns
    !That's why we have 2 instead of 1
    call RT%cmaxvol(per2, swapsmade2, maxswaps, URT)
    !Exit if no swaps were made
    if (swapsmade1 + swapsmade2 == 0) then
      exit
    end if
  end do
  !Rows should coincide with the rows of A, so we use peri
  C = Acols(Afun, M, rank, peri, per2, param)
  UR = .T.URT
  !We swap the columns back to make them coincide with the columns of A
  call UR%permcols(per2, 2)
  !Now our low-rank approximation is C \hat A^{-1} R = C*UR
end

!Use to increase number of rows and columns after maxvol
subroutine maxvol2(Afun, k, l, per1, per2, param, C, UR)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: k, l !Increased submatrix sizes
  Type(IntVec) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx) :: C !Columns for CUR approximation
  Type(Mtrx) :: UR !UR of CUR approximation

  Type(IntVec) peri !Identity permutation
  Integer(4) M, N !Matrix sizes
  Integer(4) rank !Desired rank
  Type(Mtrx) R !Matrix rows
  Type(Mtrx) Ahat !Submatrix \hat A
  Type(Mtrx) U, S, V !SVD of the submatrix
  Integer(4) i
  
  !Get the unknown sizes from input
  M = C%n
  N = UR%m
  rank = C%m
  
  !Initialize identity permutation
  call peri%perm(M)
  
  !We need to start from good rows and columns, so
  !that the top left r by r matrix has maximum volume.
  !Here we use already calculated C and A^{-1} R
  !We return back the permutations
  call C%permrows(per1, 1)
  call UR%permcols(per2, 1)
  !Can be applied to the entire matrix like A%hmaxvol2 too, if necessary.
  call C%hmaxvol2(1, rank, k, per1, 0)
  !2 indicates that we swap columns (1 if rows)
  !1 indicates that our rows are already multiplied be A^-1 (0 if not)
  call UR%hmaxvol2(2, rank, l, per2, 1)
  !We now need to use more rows and columns
  C = Acols(Afun, M, l, peri, per2, param)
  R = Arows(Afun, k, N, per1, per2, param)
  !Moreover, we need to use PROJECTIVE VOLUME
  Ahat = R%subarray(k,l)
  call Ahat%svd(U, S, V)
  !We use r-pseudoinverse
  U = U%subarray(k,rank)
  S = S%subarray(rank,rank)
  V = V%subarray(rank,l)
  do i = 1, rank
    S%d(i,i) = 1.0d0/S%d(i,i)
  end do
  call R%permcols(per2, 2)
  !Multiplication in stable order
  UR = ((.T.V)*S) * ((.T.U)*R)
end

!CUR with arbitrary number of rows and columns
!Can be used to improve maxvol/maxvol2 or standalone
subroutine maxvolproj(Afun, M, N, rank, k, l, per1, per2, param, C, UR, maxsteps, maxswaps)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: M, N !Sizes of A
  Integer(4), intent(in) :: rank !Desired rank
  Integer(4), intent(in) :: k, l !Desired number of rows and columns for CUR
  Type(IntVec) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx), intent(out) :: C !Columns for CUR approximation
  Type(Mtrx), intent(out) :: UR !UR of CUR approximation
  Integer(4), intent(in) :: maxsteps, maxswaps !Maximum number of steps and swaps

  Type(IntVec) peri !Dummy permutation
  Type(Mtrx) R !Matrix rows
  Type(Mtrx) Ahat !Submatrix \hat A
  Type(Mtrx) U, S, V !SVD of the submatrix
  Integer(4) swapsmade1, swapsmade2 !Number of swaps in rows and columns
  Integer(4) i

  !We do essentially the same we have been doing with 'cmaxvol',
  !but now we use 'dominantc' and 'dominantr'

  !maxvol-rect of size k x r: find k good rows
  do i = 1, maxsteps
    R = Arows(Afun, k, N, per1, per2, param)
    call R%dominantr(2, rank, k, per2, swapsmade1, maxswaps)
    C = Acols(Afun, M, rank, per1, per2, param)
    !1 for swaps of rows; 0 for no rows explicitly kept unswapped
    !(first rows can be saved to preserve the r x r dominant submatrix)
    call C%dominantc(1, rank, k, 0, per1, per2, swapsmade2, maxswaps)
    if (swapsmade1 + swapsmade2 == 0) then
      exit
    end if
  end do
  !We save per1 and work in the copy
  call peri%copy(per1)
  !maxvol-rect of size r x l: find l good columns
  do i = 1, maxsteps
    C = Acols(Afun, M, l, peri, per2, param)
    call C%dominantr(1, rank, l, peri, swapsmade1, maxswaps)
    R = Arows(Afun, rank, N, peri, per2, param)
    call R%dominantc(2, rank, l, 0, peri, per2, swapsmade2, maxswaps)
    if (swapsmade1 + swapsmade2 == 0) then
      exit
    end if
  end do
  !Return peri to identity permutation
  call peri%deinit()
  call peri%perm(M)
  !Projective volume business like in maxvol2
  C = Acols(Afun, M, l, peri, per2, param)
  R = Arows(Afun, k, N, per1, per2, param)
  Ahat = R%subarray(k,l)
  call Ahat%svd(U, S, V)
  U = U%subarray(k,rank)
  S = S%subarray(rank,rank)
  V = V%subarray(rank,l)
  do i = 1, rank
    S%d(i,i) = 1.0d0/S%d(i,i)
  end do
  call R%permcols(per2, 2)
  UR = ((.T.V)*S) * ((.T.U)*R)
end

!Truncated SVD on CUR
!Decreases approximation rank (therefore, use higher rank in advance)
!ALWAYS use for geometrically (exponentially) decreasing singular values!!!
!(allows to get close to exact SVD much much faster)
subroutine TruncateCUR(C, UR, new_rank, err_bound)
  Type(Mtrx) :: C, UR !Factors of CUR approximation
  Integer(4), intent(in) :: new_rank !Maximum allowed rank
  Double precision, intent(in), optional :: err_bound !Error bound
  !Nullifies all singular values smaller than err_bound
  
  Double precision tau !Another name for err_bound
  Integer(4) rank !Another name for new_rank
  Integer(4) k !Initial rank
  Type(Mtrx) R1, L2 !Upper triangular and lower triangular matrices from QR and LQ
  Type(Mtrx) Q1, Q2 !Q-factors for QR and LQ
  Type(Vector) tau1, tau2 !Information about Q1 and Q2
  Type(Mtrx) U, S, V !For SVD
  Type(Mtrx) M !Just a temporary matrix
  Integer(4) i

  k = C%m
  
  if (present(err_bound)) then
    tau = err_bound
  else
    tau = 0.0d0
  end if
  !Make rank of appropriate size if necessary
  rank = max(1, new_rank)
  rank = min(rank, k)
  
  call C%halfqr(Q1, tau1, R1)
  call UR%halflq(L2, tau2, Q2)
  
  M = R1*L2
  
  call M%svd(U, S, V)
  do i = 1, rank
    if (S%d(i,i) <= tau) then
      rank = i-1
      exit
    end if
  end do
  
  !Exit if no truncation needed
  if (rank == k) then
    return
  end if
  
  !Multiply Q1 and Q2 by svd factors to get new approximation
  U = U%subarray(k, rank)
  V = S%subarray(rank, rank) * V%subarray(rank, k)
  C = U%multq(Q1, tau1, 'L', 'D')
  UR = V%multq(Q2, tau2, 'R', 'U')
end

end module
