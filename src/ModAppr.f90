Module ModAppr
  USE ModSparse
  
  !Module for constructing fast approximations
  !Contains adaptive cross (ACA), maxvol, maxvol2, maxvolproj and TruncateCUR

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
    !$OMP PARALLEL DO
    do i = 1, N
      res%d(i,j) = Afun(per1%d(i), per2%d(j), param)
    end do
    !$OMP END PARALLEL DO
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
    !$OMP PARALLEL DO
    do i = 1, N
      res%d(j,i) = Afun(per1%d(i), per2%d(j), param)
    end do
    !$OMP END PARALLEL DO
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
    !$OMP PARALLEL DO
    do j = 1, N
      res%d(i,j) = Afun(per1%d(i), per2%d(j), param)
    end do
    !$OMP END PARALLEL DO
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
    !$OMP PARALLEL DO
    do i = 1, N
      res%d(i,j) = Afun(per1%d(j), per2%d(i), param)
    end do
    !$OMP END PARALLEL DO
  end do
end

!Adaptive Cross. Improved code of Stanislav Stavtsev. Unlike versions by Stavtsev, Savostianov and Bebendorf,
!uses rook pivoting instead of column pivoting, extends it to more starting columns and adaptive starting columns
!and allows rho-locally maximum search.
subroutine ACA(Afun, param, Ni, Nj, MaxRank, jpmax_, U, V, per, rho_, rel_err, abs_err)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Type(Mtrx), intent(in) :: param !Matrix parameters
  Integer(4), intent(in) :: Ni, Nj !Matrix sizes
  Integer(4), intent(in) :: MaxRank !Maximum rank
  Integer(4), intent(in) :: jpmax_ !Number of starting columns; 0 is adaptive
  Type(Mtrx), intent(out) :: U, V !Output low-rank factors
  Type(IntVec), intent(in) :: per !Permutation to select starting rows (should be random)
  Double precision, intent(in), optional :: rho_ !rho-locally maximum
  Double precision, intent(in), optional :: rel_err !Desired relative error
  Double precision, intent(in), optional :: abs_err !Desired absolute error; Max(rel,abs) is chosen
       
  Type(Mtrx) Up !Starting columns
  Integer(4) jpmax !jpmax used
  Logical adaptive !Whether number of starting rows is increased each iteration
       
  Logical, allocatable :: KnC(:) !Used columns toggle
  Integer(4), allocatable :: KpC(:) !Starting columns toggle with indices
  Integer(4), allocatable :: jPs(:) !Starting columns indices
  Integer(4), allocatable :: iNs(:), jNs(:) !Used rows and columns indices
  Integer(4) jP !New starting column index
  Integer(4) iN, jN !New row and column indices
  Integer(4) NSkel !Current skeletion (submatrix) size
  Double precision Elem, ElemS !Pivot element and its square root
  Double precision MaxElC, MaxElR !Max element in column and row
  Double precision rho !rho used
  Double precision err_bound !Error bound
  
  !Double precision, allocatable :: uvec(:), vvec(:) !Vectors for F-norm update
  !Double precision fnorm !Current Frobenius norm
  !Currently not used, because too expensive
       
  Double precision, allocatable :: tmp_array(:,:) !For copying
  
  Integer(4) i, j, s !Indices
       
  if ((Ni <= 0) .or. (Nj <= 0) .or. (maxrank <= 0) .or. (jpmax_ < 0)) then
    print *, 'Input sizes to ACA should be positive!'
    return
  end if
       
  rho = 1000000000.0d0
  if (present(rho_)) then
    if (rho_ < 1) then
      print *, 'Input rho in ACA should be at least 1!'
      return
    end if
    rho = rho_
  end if
       
  if (present(abs_err)) then
    err_bound = abs_err
  else
    err_bound = 0
  end if

  jpmax = jpmax_
  if (jpmax_ == 0) then
    jpmax = 1
    adaptive = .true.
  else
    adaptive = .false.
  end if
  
  if ((adaptive) .and. (maxrank*2 > Nj)) then
    print *, 'Adaptive starting columns number in ACA may exceed total number of columns!'
  end if
  if (maxrank + jpmax > Nj) then
    print *, 'Number of used and starting columns in ACA may exceed total number of columns!'
  end if
  if (maxrank > min(Ni,Nj)) then
    print *, 'Maximum rank in ACA exceeds total number of rows or columns!'
  end if

  call U%init(Ni,maxrank)
  call V%init(Nj,maxrank)
  if (adaptive) then
    call Up%init(Ni,maxrank)
    Up%m = 1
    Allocate(jPs(maxrank))
  else
    Allocate(jPs(jpmax))
    call Up%init(Ni,jpmax)
  end if
  Allocate(KnC(Nj))
  Allocate(KpC(Nj))
  Allocate(iNs(maxrank))
  Allocate(jNs(maxrank))

  KnC = .false.
  KpC = 0
  jPs = 0
  iNs = 0
  jNs = 0
        
  !Init column indices 
  jN = 1
        
  !Init row index
  iN = 1
        
  do j = 1, jpmax
    KpC(per%d(j)) = j
    jPs(j) = per%d(j)
  end do
  jP = jpmax+1
        
  NSkel = 0
        
  !Initialize starting columns
  do j = 1, jpmax
    do i = 1, Ni
      Up%d(i,j) = Afun(i, per%d(j), param)
    end do
  end do

  !Start approximation process
  do while (NSkel < MaxRank)
        
    !Number of skeletons to aproximate matrix
    NSkel = NSkel + 1    
            
    !Choose new starting column if it was used
    if ((KpC(jN) > 0) .and. (Nskel > 1)) then
      s = KpC(jN)
      KpC(jN) = 0
      do while (KnC(per%d(jP)))
        jP = jP + 1
      end do
      KpC(per%d(jP)) = s
      jPs(s) = per%d(jP)
      do i = 1, Ni
        Up%d(i,s) = Afun(i,per%d(jP),param)
      end do
      do j = 1, NSkel-2
        Up%d(:,s) = Up%d(:,s) - U%d(:,j)*V%d(per%d(jP),j)
      end do
      Up%d(iNs(1:Nskel-2),s) = 0
      jP = jP + 1
    end if
            
    !Adaptive increase of the number of starting columns
    if ((adaptive) .and. (Nskel > 1)) then
      jpmax = jpmax + 1
      s = jpmax
      Up%m = jpmax
      do while (KnC(per%d(jP)))
        jP = jP + 1
      end do
      KpC(per%d(jP)) = s
      jPs(s) = per%d(jP)
      do i = 1, Ni
        Up%d(i,s) = Afun(i,per%d(jP),param)
      end do
      do j = 1, NSkel-2
        Up%d(:,s) = Up%d(:,s) - U%d(:,j)*V%d(per%d(jP),j)
      end do
      Up%d(iNs(1:Nskel-2),s) = 0
      jP = jP + 1
    end if
                        
    !Update elements in starting columns
    if (NSkel > 1) then
      do j = 1, jpmax
        Up%d(:,j) = Up%d(:,j) - U%d(:,NSkel-1) * V%d(jPs(j),NSkel-1)
        Up%d(iN,j) = 0
      end do
    end if

    !Find the maximum element in starting columns
    call Up%cnormloc(s,jN)
            
    !Loop until rho-maximum is reached
    MaxElC = Up%d(s,jN)
    MaxElR = 0
    do while (abs(MaxElC) > rho*abs(MaxElR))
            
      iN = s
            
      !Calculate elements in new row
      do j = 1, Nj
        V%d(j,NSkel) = Afun(iN,j,param)
      end do
      do j = 1, NSkel-1
        V%d(:,NSkel) = V%d(:,NSkel) - U%d(iN,j)*V%d(:,j)
      end do
      V%d(jNs(1:Nskel-1),NSkel) = 0

      !Calculate the maximum in the new row
      s = jN
      jN = myidamax(Nj, V%d(:,Nskel))
      MaxElR = V%d(jN,NSkel)
            
      !Do not calculate the column, if it's the same
      if ((s .ne. jN) .or. (U%d(iN,NSkel) == 0)) then
            
        !Calculate elements in new column
        if (KpC(jN) == 0) then
          do i = 1, Ni
            U%d(i,NSkel) = Afun(i,jN,param)
          end do
          do j = 1, NSkel-1
            U%d(:,NSkel) = U%d(:,NSkel) - U%d(:,j)*V%d(jN,j)
          end do
          U%d(iNs(1:Nskel-1),NSkel) = 0
        else
          U%d(:,NSkel) = Up%d(:,KpC(jN))
        end if
            
        !Calculate the maximum in the new column
        s = myidamax(Ni, U%d(:,Nskel))
        MaxElC = U%d(s,NSkel)
      end if
              
    end do
            
    !Check error bound, exit if reached
    if (present(rel_err)) then
      err_bound = max(err_bound, rel_err*abs(MaxElC))
    end if
    if (NSkel > 1) then
      if (abs(MaxElC) < err_bound) then
        Nskel = NSkel-1
        exit
      end if
    end if
            
    !Form skeleton at cross iN, jN
    ElemS = sqrt(abs(MaxElR))
    Elem = Elems/MaxElR
    ElemS = 1/ElemS
    if (MaxElR .ne. 0) then
      U%d(:,NSkel) = Elem*U%d(:,NSkel)
      V%d(:,NSkel) = ElemS*V%d(:,NSkel)
    else
      print *, 'Zero Skeleton encountered!'
      NSkel = NSkel-1
      exit
    end if 

    !Mark new row and column
    KnC(jN) = .true.
    iNs(NSkel) = iN
    jNs(NSkel) = jN

  end do
  
  if (NSkel < maxrank) then
    Allocate(tmp_array(Ni,NSkel))
    call dlacpy('A', Ni, NSkel, U%d, Ni, tmp_array, Ni)
    U = tmp_array
    Allocate(tmp_array(Nj,NSkel))
    call dlacpy('A', Nj, NSkel, V%d, Nj, tmp_array, Nj)
    V = tmp_array
  end if
  
  Deallocate(KnC,KpC,jPs,iNs,jNs)
end
  
!Simplest CUR approximation with U = \hat A^{-1}
subroutine maxvol(Afun, M, N, rank, per1, per2, param, C, UR, maxsteps, maxswaps, CA, pre)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: M, N, rank !Sizes of A and the desired rank
  Type(IntVec) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx), intent(out) :: C !Columns for CUR approximation
  Type(Mtrx), intent(out) :: UR !UR of CUR approximation
  Integer(4), intent(in) :: maxsteps, maxswaps !Maximum number of steps and swaps
  Type(Mtrx), intent(out), optional :: CA !C*Ahat^{-1} or C*Rhat^{-1} from Ahat=Qhat*Rhat, whichever is faster
  Logical, intent(in), optional :: pre !Run premaxvol? [FALSE]
  Logical :: pre_

  Type(Mtrx) ABout
  Type(Mtrx) Ahat, Q, R
  Type(Vector) tau
  Type(IntVec) peri !Identity permutation
  Type(Mtrx) RT !Transposed rows
  Type(Mtrx) URT !Transposed UR
  Integer(4) swapsmade1, swapsmade2 !Number of swaps in rows and columns
  Integer(4) i
  
  if (present(pre)) then
    pre_ = pre
  else
    pre_ = .false.
  end if
  
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
    if ((i == 1) .and. (pre_)) then
      C = Acolst(Afun, M, rank, per1, per2, param)
      call C%premaxvol(rank, per1)
      !C = .T.C
      !if (present(CA)) then
      !  call C%cmaxvol(per1, swapsmade1, maxswaps, CA, ABin = ABout)
      !else
      !  call C%cmaxvol(per1, swapsmade1, maxswaps, ABin = ABout)
      !end if
      swapsmade1 = rank
    else
      C = Acols(Afun, M, rank, per1, per2, param)
      !We use column version of maxvol.
      if (present(CA)) then
        call C%cmaxvol(per1, swapsmade1, maxswaps, CA)
      else
        call C%cmaxvol(per1, swapsmade1, maxswaps)
      end if
    end if
    !Exit if no swaps were made
    if ((swapsmade1 == 0) .and. (i > 1)) then
      exit
    end if
    
     if ((i == 1) .and. (pre_)) then
       RT = Arows(Afun, rank, N, per1, per2, param)
       call RT%premaxvol(rank, per2, ABout)
       RT = .T.RT
       call RT%cmaxvol(per2, swapsmade2, maxswaps, URT, ABin = .T.ABout)
       swapsmade2 = swapsmade2 + rank
     else
      !Select first r rows (transposed)
      RT = Arowst(Afun, rank, N, per1, per2, param)
      !We again use column version and swap columns
      !That's why we have 2 instead of 1
      call RT%cmaxvol(per2, swapsmade2, maxswaps, URT)
     end if
    !Exit if no swaps were made
    if ((swapsmade2 == 0) .and. ((i > 1) .or. (.not.(pre_)))) then
      exit
    end if
  end do
  !Rows should coincide with the rows of A, so we use peri
  C = Acols(Afun, M, rank, peri, per2, param)
  if (present(CA) .and. ((swapsmade2 > 0) .or. ((i == 1) .and. (pre_)))) then
    CA = Acols(Afun, M, rank, per1, per2, param)
    Ahat = CA%subarray(rank, rank)
    call Ahat%halfqr(Q, tau, R)
    call dtrsm('R', 'U', 'N', 'N', CA%n, CA%m, 1.0d0, R%d, R%n, CA%d, CA%n)
  end if
  UR = .T.URT
  !We swap the columns back to make them coincide with the columns of A
  call UR%permcols(per2, 2)
  !Now our low-rank approximation is C \hat A^{-1} R = C*UR
end

!Use to increase number of rows and columns after maxvol
subroutine maxvol2(Afun, k, l, per1, per2, param, C, UR, ca)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Integer(4), intent(in) :: k, l !Increased submatrix sizes
  Type(IntVec) :: per1, per2 !Permutations of rows and columns
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Type(Mtrx) :: C !Columns for CUR approximation
  Type(Mtrx) :: UR !UR of CUR approximation. Not permuted!
  Logical, intent(in), optional :: ca !C is actually C Ahat^-1 (or, at least, C Rhat^-1, Ahat = Qhat Rhat)? Permuted!

  Logical :: ca_
  Type(IntVec) peri !Identity permutation
  Integer(4) M, N !Matrix sizes
  Integer(4) rank !Desired rank
  Type(Mtrx) R !Matrix rows
  Type(Mtrx) Ahat !Submatrix \hat A
  Type(Mtrx) U, V !SVD of the submatrix
  Type(Vector) s !Singular values
  
  if (present(ca)) then
    ca_ = ca
  else
    ca_ = .false.
  end if
  
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
  if (.not. ca_) then
    call C%permrows(per1, 1)
  end if
  
  call UR%permcols(per2, 1)
  !Can be applied to the entire matrix like A%umaxvol2 too, if necessary.
  call C%umaxvol2(1, rank, k, per1, ca_)
  !2 indicates that we swap columns (1 if rows)
  !.true. indicates that our rows are already multiplied by A^-1 (.false. if not)
  call UR%umaxvol2(2, rank, l, per2, .true.)
  !We now need to use more rows and columns
  C = Acols(Afun, M, l, peri, per2, param)
  R = Arows(Afun, k, N, per1, per2, param)
  !Moreover, we need to use PROJECTIVE VOLUME
  Ahat = R%subarray(k,l)
  call Ahat%svd(U, s, V)
  !We use r-pseudoinverse
  U = U%subarray(k,rank)
  s = s%subarray(rank)
  V = V%subarray(rank,l)
  call R%permcols(per2, 2)
  C = C .dT. V
  UR = (U .dd. s) .Td. R
  !Multiplication in stable order
  !UR = ((.T.V)*S) * ((.T.U)*R)
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
  Type(Mtrx) U, V !SVD of the submatrix
  Type(Vector) s !Singular values
  Integer(4) swapsmade1, swapsmade2 !Number of swaps in rows and columns
  Integer(4) i, j
  
  Type(Mtrx) ABout, CMout, Zout
  Type(Vector) cout, Lout
  
  !Type(Mtrx) Q1
  !Type(Vector) tau1

  !We do essentially the same we have been doing with 'cmaxvol',
  !but now we use 'dominantc' and 'dominantr'

  !maxvol-rect of size k x r: find k good rows
  do i = 1, maxsteps
    if (k >= l) then
      R = Arows(Afun, k, N, per1, per2, param)
      if (i == 1) then
        call R%premaxvol(rank, per2, ABout, cout)
        call R%dominantr(2, rank, k, per2, swapsmade1, maxswaps, ABin = ABout, cin = cout)
        swapsmade1 = swapsmade1 + rank
      else
        call R%dominantr(2, rank, k, per2, swapsmade1, maxswaps)
      end if
    else
      C = Acolst(Afun, M, l, per1, per2, param)
      if (i == 1) then
        call C%premaxvol(rank, per1, ABout, cout)
        call C%dominantr(2, rank, l, per1, swapsmade1, maxswaps, ABin = ABout, cin = cout)
        swapsmade1 = swapsmade1 + rank
      else
        call C%dominantr(2, rank, l, per1, swapsmade1, maxswaps)
      end if
    end if
    if ((swapsmade1 == 0) .and. (i > 1)) then
      exit
    end if
    if (k >= l) then
     if (i == 1) then
       C = Acolst(Afun, M, rank, per1, per2, param)
       call C%premaxvol(rank, per1, ABout)
       ABout%d(1:rank,1:rank) = 0
       do j = 1, rank
         About%d(j,j) = 1.0d0
       end do
       call ABout%umaxvol2(2, rank, k, per1, .true., C, CMout, Zout, Lout)
       call C%dominantc(2, rank, k, per1, swapsmade2, maxswaps, CMout, Zout, Lout)
       swapsmade2 = swapsmade2 + k
     else
      C = Acols(Afun, M, rank, per1, per2, param)
      !1 for swaps of rows; 0 for no rows explicitly kept unswapped
      !(first rows can be saved to preserve the r x r dominant submatrix)
      call C%dominantc(1, rank, k, per1, swapsmade2, maxswaps)
     end if
    else
      if (i == 1) then
       R = Arows(Afun, rank, N, per1, per2, param)
       call R%premaxvol(rank, per2, ABout)
       ABout%d(1:100,1:100) = 0
       do j = 1, 100
         About%d(j,j) = 1.0d0
       end do
       call ABout%umaxvol2(2, rank, l, per2, .true., R, CMout, Zout, Lout)
       call R%dominantc(2, rank, l, per2, swapsmade2, maxswaps, CMout, Zout, Lout)
       swapsmade2 = swapsmade2 + l
     else
      R = Arows(Afun, rank, N, per1, per2, param)
      call R%dominantc(2, rank, l, per2, swapsmade2, maxswaps)
     end if
    end if
    if (swapsmade2 == 0) then
      exit
    end if
  end do
  
  !Instead of the second cycle, just work in Ahat^+ R: it is faster and gives better accuracy
  if (k >= l) then
    R = Arows(Afun, k, N, per1, per2, param)
    Ahat = R%subarray(k,rank)
    R = Ahat .Id. R
    call R%dominantc(2, rank, l, per2, swapsmade2, maxswaps)
  else
    C = Acols(Afun, M, l, per1, per2, param)
    Ahat = C%subarray(rank,l)
    C = C .dI. Ahat
    call C%dominantc(1, rank, k, per1, swapsmade2, maxswaps)
  end if
  
!Inverse through SVD of columns C: improve is not observable (need column diagram to see the difference)
!   call peri%perm(N)
!   C = Acols(Afun, M, l, per1, per2, param)
!   R = Arows(Afun, k, N, per1, peri, param)
!   call C%halfqr(Q1, tau1, Ahat)
!   call Ahat%svdr(rank, U, V)
!   C = U%multq(Q1, tau1, 'L', 'D')
!   Ahat = C%subarray(k, rank)
!   call C%permrows(per1, 2)
!   UR = Ahat .Id. R
  
  !Projective volume business like in maxvol2
  call peri%perm(M)
  C = Acols(Afun, M, l, peri, per2, param)
  R = Arows(Afun, k, N, per1, per2, param)
  Ahat = R%subarray(k,l)
  call Ahat%svd(U, s, V)
  U = U%subarray(k,rank)
  s = s%subarray(rank)
  V = V%subarray(rank,l)
  call R%permcols(per2, 2)
  C = C .dT. V
  UR = (U .dd. s) .Td. R
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
  Type(Mtrx) U, V !For SVD
  Type(Vector) s !Singular values
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
  
  call M%svd(U, s, V)
  do i = 1, rank
    if (s%d(i) <= tau) then
      rank = i-1
      exit
    end if
  end do
  if (present(trunc_err)) then
    trunc_err = norm2(s%d(rank+1:k))
  end if
  
  !Exit if no truncation needed
  if (rank == k) then
    return
  end if
  
  !Multiply Q1 and Q2 by svd factors to get new approximation
  U = U%subarray(k, rank)
  V = s%subarray(rank) .dot. V%subarray(rank, k)
  C = U%multq(Q1, tau1, 'L', 'D')
  UR = V%multq(Q2, tau2, 'R', 'U')
end

!Alternating projections for nonnegative matrix approximation
!with truncated SVD replaced by CUR
subroutine PositCUR(Afun, param, M, N, rank, k2, k3, minelem, maxelem, epsmult, C, AR, verbose)
  procedure(elem_fun) :: Afun !Function, returning elements of A
  Type(Mtrx), intent(in) :: param !Parameters for Afun
  Integer(4), intent(in) :: M, N, rank !Sizes of A and the desired rank
  Integer(4) :: k2, k3 !Submatrix sizes for TSVD and projective volume
  Double precision, intent(in) :: minelem !Desired minimum matrix element
  Double precision, intent(in), optional :: maxelem !Desired minimum matrix element
  Double precision, intent(in) :: epsmult !Error multiplicator for shift calculation
  Type(Mtrx), intent(out) :: C, AR !Output aprroximation of rank "rank"
  Integer(4), intent(in) :: verbose !0-2; how much to print.
  
  Type(Mtrx) C1 !To speed up maxvol2, we save C \hat A^{-1} in C, so C1 is dummy for C
  Integer(4) maxsteps, maxswaps
  Type(IntVec) peri, peri2, per1, per2
  Type(Mtrx) R1, L2, R, U, V, Ahat, Q1, Q2, param_uv, EU, EV
  Type(Vector) s
  Integer(4) i, j, steps, cur, bads
  Double precision eps2, badstot, maxelem_, elem, epsmin, tmp
  Type(Vector) tau, tau2
  
  Double precision normf
  
  if (present(maxelem)) then
    maxelem_ = maxelem
  else
    maxelem_ = minelem
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
  
  call maxvol(Afun, M, N, k2, per1, per2, param, C1, AR, maxsteps, maxswaps, C)
  call maxvol2(Afun, k3, k3, per1, per2, param, C, AR, .true.)
  
  C = Acols(Afun, M, k3, peri, per2, param)
  R = Arows(Afun, k3, N, per1, per2, param)
  Ahat = R%subarray(k3,k3)
  call Ahat%svd(U, s, V)
  cur = k2
  do i = k3, 1, -1
    if (((s%d(i) > eps*k3*s%d(1)) .or. (i <= rank)) .and. (i <= k2)) then
      s%d(i) = 1.0d0/s%d(i)
    else
      s%d(i) = 0
      cur = i-1
    end if
  end do
  call R%permcols(per2, 2)
  U = U%subarray(k3,cur)
  V = s%subarray(cur) .dot. V%subarray(cur,k3)
  C = C .dT. V
  AR = U .Td. R
  
  if (cur > rank) then
    call C%halfqr(Q1, tau, R1)
    call AR%halflq(L2, tau2, Q2)
    Ahat = R1*L2
    call Ahat%svd(U, s, V)
    U = U%subarray(cur, rank)
    V = s%subarray(rank) .dot. V%subarray(rank, cur)
    C = U%multq(Q1, tau, 'L', 'D')
    AR = V%multq(Q2, tau2, 'R', 'U')
  else
    call C%halfqr(Q1, tau, R1)
    call AR%halflq(L2, tau2, Q2)
    Ahat = R1*L2
    call Ahat%svd(U, s, V)
  end if
  normf = norm2(s%d(1:rank))
  
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
  param_uv%d(2,2*rank+1) = maxelem_
  
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
    call Ahat%svd(U, s, V)
    !Replace sqrt(eps) with eps? It got bad errors sometimes, but can't find why and when.
    do i = k3, 1, -1
      if (((s%d(i) > eps*k3*s%d(1)) .or. (i <= rank+1)) .and. (i <= k2)) then
        s%d(i) = 1.0d0/s%d(i)
      else
        s%d(i) = 0
        cur = i
      end if
    end do
    call R%permcols(per2, 2)
    U = U%subarray(k3,cur)
    V = s%subarray(cur) .dot. V%subarray(cur,k3)
    C = C .dT. V
    AR = U .Td. R
  
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
    param_uv%d(2,2*rank+1) = maxelem_-eps2
  
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
      R = EV .dT. R1
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
        if ((elem < minelem) .or. ((elem > maxelem_) .and. (maxelem_ .ne. minelem))) then
          bads = bads + 1
          exit
        end if
      end do
      if (bads > 0) then
        exit
      end if
    end do
    if ((bads == 0) .or. (steps > 200)) then
      Ahat = Arows(Afun, M, N, peri, peri2, param)
      U = Ahat - C*AR
      U = U*(1.0d0/Ahat%fnorm())
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
