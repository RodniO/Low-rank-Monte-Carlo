Module ModAppr
  USE ModSparse
  
  !Module for constructing fast approximations
  !Contains maxvol, maxvol2, maxvolproj and TruncateCUR

  abstract interface
    !Interface for function, returning matrix elements
    pure function elem_fun(i, j, param) Result(res)
      USE ModMtrx
      Integer(4), intent(in) :: i, j
      Type(Mtrx), intent(in) :: param
      Double precision :: res
    end function
  end interface
  
  Contains
  
!Calculates matrix element for i-th row and j-th column.
pure function Afun_elem(i, j, param) Result(res)
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  !Just uses parameters array as the input matrix
  res = param%d(i,j)
end

pure function Afun_uv(i, j, param) Result(res)
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  
  Integer(4) k, i1
  
  k = param%m/2
  res = 0.0d0
  do i1 = 1, k
    res = res + param%d(i,i1)*param%d(j,i1+k)
  end do
end

pure function Afun_uv_err(i, j, param) Result(res)
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  
  Integer(4) k, i1
  
  k = param%m/3
  res = 0.0d0
  do i1 = 1, k
    res = res + param%d(i,i1)*param%d(j,i1+k)
  end do
end

pure function Afun_uvmod(i, j, param) Result(res)
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  
  Integer(4) k, i1
  
  k = param%m/2
  res = 0.0d0
  do i1 = 1, k
    res = res + param%d(i,i1)*param%d(j,i1+k)
  end do
  if (res < param%d(1,2*k+1)) then
    res = param%d(1,2*k+1)
  else if ((res > param%d(2,2*k+1)) .and. (param%d(2,2*k+1) > param%d(1,2*k+1))) then
    res = param%d(2,2*k+1)
  end if
end

pure function Afun_uvmod_err(i, j, param) Result(res)
  Integer(4), intent(in) :: i, j !row and column indices
  Type(Mtrx), intent(in) :: param !Arbitrary array of parameters
  Double precision :: res !Output value
  
  Integer(4) k, i1
  Double precision res2, resmin, resmax
  
  k = param%m/3
  res = 0.0d0
  do i1 = 1, k
    res = res + param%d(i,i1)*param%d(j,i1+k)
  end do
  res2 = 0.0d0
  do i1 = 2*k+1, 2*k+k/2
    res2 = res2 + param%d(i,i1)*param%d(j,i1+k/2)
  end do
  resmin = param%d(1,3*k+1) + max(0.0d0, res2)
  resmax = param%d(2,3*k+1) + min(0.0d0, res2)
  if (res < resmin) then
    res = resmin
  else if ((res > resmax) .and. (param%d(2,3*k+1) > param%d(1,3*k+1) )) then
    res = resmax
  end if
end
  
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

!Returns r columns of length N with elements from Afun
function Acolst(Afun, N, r, per1, per2, param) Result(res)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: N, r !Length and number of columns
  Type(IntVec), intent(in) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx) :: res !Output matrix of columns
  
  Integer(4) i, j
  
  call res%init(r, N)
  do j = 1, r
    do i = 1, N
      res%d(j,i) = Afun(per1%d(i), per2%d(j), param)
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
subroutine TruncateCUR(C, UR, new_rank, err_bound, trunc_err)
  Type(Mtrx) :: C, UR !Factors of CUR approximation
  Integer(4), intent(in) :: new_rank !Maximum allowed rank
  Double precision, intent(in), optional :: err_bound !Error bound
  Double precision, intent(out), optional :: trunc_err !Truncation error
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
  if (present(trunc_err)) then
    trunc_err = 0
    do i = rank+1, k
      trunc_err = trunc_err + S%d(i,i)**2
    end do
    trunc_err = sqrt(trunc_err)
  end if
  
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

!Alternating projections for nonnegative matrix approximation
!with truncated SVD replaced by CUR
subroutine PositCUR(Afun, param, M, N, rank, k2, k3, minelem, maxelemin, epsmult, C, AR, verbose)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Integer(4), intent(in) :: M, N, rank !Sizes of A and the desired rank
  Integer(4) :: k2, k3 !Submatrix sizes for TSVD and projective volume
  Double precision, intent(in) :: minelem !Desired minimum matrix element
  Double precision, intent(in), optional :: maxelemin !Desired minimum matrix element
  Double precision, intent(in) :: epsmult !Error multiplicator for shift calculation
  Type(Mtrx), intent(out) :: C, AR !Output aprroximation of rank "rank"
  Integer(4), intent(in) :: verbose !0-2; how much to print.
  
  Integer(4) maxsteps, maxswaps
  Type(IntVec) peri, peri2, per1, per2
  Type(Mtrx) R1, L2, R, U, S, V, Ahat, Q1, Q2, param_uv, EU, EV
  Integer(4) i, j, steps, cur, bads
  Double precision eps2, badstot, maxelem, elem, epsmin, tmp
  Type(Vector) tau, tau2
  
  Double precision normf
  
  if (present(maxelemin)) then
    maxelem = maxelemin
  else
    maxelem = minelem
  end if
  
  !Maximum number of steps for maxvol
  maxsteps = 4
  !Maximum number of row and column swaps for maxvol
  maxswaps = 2*rank

  call per1%perm(M)
  call per2%perm(N)
  call peri%perm(M)
  call peri2%perm(N)
  
  !Start should be good, so random rows and columns are chosen and then improved with premaxvol
  do i = 1, k3
    call random_number(tmp)
    j = k3+1 + FLOOR((M-k3)*tmp)
    call per1%swap(i,j)
  end do
  R = Arows(Afun, k3, N, per1, per2, param)
  call R%premaxvol(k2, per2)
  Ahat = .T.(R%subarray(k3,k2))
  call Ahat%premaxvol(k2, per1)
  
  call maxvol(Afun, M, N, k2, per1, per2, param, C, AR, maxsteps, maxswaps)
  call maxvol2(Afun, k3, k3, per1, per2, param, C, AR)
  
  C = Acols(Afun, M, k3, peri, per2, param)
  R = Arows(Afun, k3, N, per1, per2, param)
  Ahat = R%subarray(k3,k3)
  call Ahat%svd(U, S, V)
  cur = k2
  do i = k3, 1, -1
    if (((S%d(i,i) > eps*k3*S%d(1,1)) .or. (i <= rank)) .and. (i <= k2)) then
      S%d(i,i) = 1.0d0/S%d(i,i)
    else
      S%d(i,i) = 0
      cur = i-1
    end if
  end do
  call R%permcols(per2, 2)
  U = U%subarray(k3,cur)
  V = S%subarray(cur,cur)*V%subarray(cur,k3)
  C = C*(.T.V)
  AR = (.T.U)*R
  
  if (cur > rank) then
    call C%halfqr(Q1, tau, R1)
    call AR%halflq(L2, tau2, Q2)
    Ahat = R1*L2
    call Ahat%svd(U, S, V)
    normf = 0
    do i = 1, rank
      normf = normf + S%d(i,i)**2
    end do
    normf = sqrt(normf)
    U = U%subarray(cur, rank)
    V = S%subarray(rank, rank) * V%subarray(rank, cur)
    C = U%multq(Q1, tau, 'L', 'D')
    AR = V%multq(Q2, tau2, 'R', 'U')
  else
    call C%halfqr(Q1, tau, R1)
    call AR%halflq(L2, tau2, Q2)
    Ahat = R1*L2
    call Ahat%svd(U, S, V)
    normf = 0
    do i = 1, rank
      normf = normf + S%d(i,i)**2
    end do
    normf = sqrt(normf)
  end if
  
  call param_uv%init(max(M,N),2*rank+1)
  do i = 1, rank
    do j = 1, M
      param_uv%d(j,i) = C%d(j,i)
    end do
    do j = 1, N
      param_uv%d(j,i+rank) = AR%d(i,j)
    end do
  end do
  param_uv%d(1,2*rank+1) = minelem
  param_uv%d(2,2*rank+1) = maxelem
  
  call EU%init(M,rank*2)
  call EV%init(N,rank*2)
  do i = 1, rank
    do j = 1, M
      EU%d(j,i+rank) = C%d(j,i)
    end do
    do j = 1, N
      EV%d(j,i+rank) = -AR%d(i,j)
    end do
  end do
  epsmin = 0

  badstot = M*N
  
  if (verbose > 0) then
    print *, 'Initial approximation done.'
  end if
  
  !k2 = k3
  
  maxsteps = 2
  steps = 0
  do while (.true.)
    steps = steps + 1
    
    C = Acols(Afun_uvmod, M, k3, peri, per2, param_uv)
    R = Arows(Afun_uvmod, k3, N, per1, per2, param_uv)
    
    Ahat = R%subarray(k3,k3)
    call Ahat%svd(U, S, V)
    !Replace sqrt(eps) with eps? It got bad errors sometimes, but can't find why and when.
    do i = k3, 1, -1
      if (((S%d(i,i) > eps*k3*S%d(1,1)) .or. (i <= rank+1)) .and. (i <= k2)) then
        S%d(i,i) = 1.0d0/S%d(i,i)
      else
        S%d(i,i) = 0
        cur = i
      end if
    end do
    call R%permcols(per2, 2)
    U = U%subarray(k3,cur)
    V = S%subarray(cur,cur)*V%subarray(cur,k3)
    C = C*(.T.V)
    AR = (.T.U)*R
  
    call TruncateCUR(C, AR, rank)
    eps2 = epsmult*epsmin/sqrt(badstot)
    do i = 1, rank
      do j = 1, M
        param_uv%d(j,i) = C%d(j,i)
      end do
      do j = 1, N
        param_uv%d(j,i+rank) = AR%d(i,j)
      end do
    end do
    param_uv%d(1,2*rank+1) = minelem+eps2
    param_uv%d(2,2*rank+1) = maxelem-eps2
  
    if (epsmin == 0.0d0) then
      do i = 1, rank
        do j = 1, M
          EU%d(j,i) = C%d(j,i)
        end do
        do j = 1, N
          EV%d(j,i) = AR%d(i,j)
        end do
      end do
      call EU%halfqr(Q1, tau, R1)
      R = EV*(.T.R1)
      epsmin = max(R%fnorm(),normf*eps/10.0d0*sqrt(badstot))
    
      if (verbose > 1) then
        print *, '1) Distance to rank r 2) shift:'
        print *, epsmin, eps2
      end if
    end if
    
    if (verbose > 1) then
      print *, 'Step', steps
    end if
  
    !call maxvol(Afun_uvmod, M, N, rank, per1, per2, param_uv, C, AR, maxsteps, maxswaps)
    
    bads = 0
    do j = 1, N
      R = C*AR%subarray(rank, per2%d(j), 1, per2%d(j))
      do i = 1, M
        elem = R%d(per1%d(i),1)
        if ((elem < minelem) .or. ((elem > maxelem) .and. (maxelem .ne. minelem))) then
          bads = bads + 1
          exit
        end if
      end do
      if (bads > 0) then
        exit
      end if
    end do
    if ((bads == 0) .or. (steps > 200)) then
      S = Arows(Afun, M, N, peri, peri2, param)
      U = S - C*AR
      U = U*(1.0d0/S%fnorm())
      exit
    end if
    
    !premaxvol не должен портить maxvol!!!
    do i = rank+1, k3
      call random_number(tmp)
      j = rank+1 + FLOOR((M-rank)*tmp)
      call per1%swap(i,j)
    end do
    R = Arows(Afun_uvmod, k3, N, per1, per2, param_uv)
    call R%premaxvol(k2, per2)
    do i = k2+1, k3
      call random_number(tmp)
      j = k2+1 + FLOOR((N-k2)*tmp)
      call per2%swap(i,j)
    end do
    C = Acolst(Afun_uvmod, M, k3, per1, per2, param_uv)
    call C%premaxvol(k2, per1)
    call maxvol(Afun_uvmod, M, N, rank, per1, per2, param_uv, C, AR, maxsteps, maxswaps)
    
  end do
  
  if (verbose == 1) then
    print *, 'STEPS:', steps
  end if
  
  if (bads > 0) then
    print *, 'DID NOT CONVERGE!'
  else
    if (verbose > 0) then
      print *, 'Nonnegative approximation finished.'
    end if
  end if
end

end module
