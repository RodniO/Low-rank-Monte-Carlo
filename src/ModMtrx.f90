Module ModMtrx
  Use ModVec

  !Wrapper for LAPACK and BLAS matrix subroutines
  
  !Also contains dominant submatrix search
  !(premaxvol, maxvol, maxvol2, dominantr, dominantc)
  ![1] A.I. Osinsky. Rectangular maximum volume and projective volume search algorithms // arXiv 1809.02334 (Submitted on 7 Sep 2018)
  
  !And some other algorithms described in
  ![2] M. Gu, S. C. Eisenstat, efficient algorithms for computing a strong rank-revealing qr factorization // SIAM J. ScI. COMPUT. — 1996. — Vol. 17, no. 4. — P. 848-869.
  ![3] How to find a good submatrix / S.A. Goreinov, I.V. Oseledets, D.V. Savostyanov et al. // Matrix Methods: Theory, Algorithms, Applications / Ed. by V. Olshevsky, E. Tyrtyshnikov. — World Scientific Publishing, 2010. — P. 247-256.
  ![4] Michalev A.Y., Oseledets I.V. Rectangular maximum-volume submatrices and their applications // Linear Algebra and its Applications. — 2018. — Vol. 538. — P. 187–211.
  ![5] A.I. Osinsky, N.L. Zamarashkin. Pseudo-skeleton approximations with better accuracy estimates // Linear Algebra and its Applications — 2018. — Vol. 537. — P. 221-249.

  !Matrix type
  Type Mtrx
    Integer(4) n !Number of rows
    Integer(4) m !Number of columns
    DOUBLE PRECISION,Allocatable :: d(:,:) !Array of matrix elements of size n by m
    Contains
      Procedure :: init => Mtrx_constructor !Allocates and initializes with zeros
      Procedure :: deinit => Mtrx_destructor !Deallocates
      Procedure :: gauss => mtrx_gauss !Random gaussian matrix
      Procedure :: random => mtrx_random !Random orthogonal matrix
      Procedure :: mask => mtrx_mask !Random mask of 0-s and 1-s
      Procedure :: badrandom => mtrx_badrandom !Random matrix with uniform elements
      Procedure :: subarray => mtrx_subarray !Returns submatrix
      Procedure :: subcol => mtrx_subcol !Returns subset of columns
      Procedure :: svd => mtrx_svd !Performs full SVD
      Procedure :: rtsolve => mtrx_rtsolve !Right triangular solver
      Procedure :: ltsolve => mtrx_ltsolve !Left triangular solver
      Procedure :: halfqr => mtrx_halfqr !Finds R of QR decomposition
      Procedure :: qr => mtrx_qr !Performs QR decomposition
      Procedure :: halflq => mtrx_halflq !Fands L of LQ decomposition
      Procedure :: lq => mtrx_lq !Performs LQ decomposition
      Procedure :: multq => mtrx_multq !Multiplies by a set of orthogonal reflectors
      Procedure :: hqr => mtrx_hqr !Householder QR with pivoting
      Procedure :: geqr => mtrx_geqr !Gu and Eisenstatt's QR [2]
      Procedure :: tauinverse => mtrx_tauinverse !tau-pseudoinverse
      Procedure :: swap => mtrx_swap !Swaps two rows or two columns
      Procedure :: permrows => mtrx_permrows !Permutes rows
      Procedure :: permcols => mtrx_permcols !Permutes columns
      Procedure :: cmaxvol => mtrx_cmaxvol, mtrx_cmaxvolro !maxvol in columns [3]
      Procedure :: fbmaxvol => mtrx_fbmaxvol !Gaussian elimination with complete pivoting
      Procedure :: maxvol2 => mtrx_maxvol2 !Original maxvol2 (aka rect_maxvol) [4]
      Procedure :: hmaxvol2 => mtrx_hmaxvol2 !Faster (Householder-based) maxvol2 [1]
      Procedure :: dominantc => mtrx_dominantc !Dominant-C [1]
      Procedure :: premaxvol => mtrx_premaxvol !Pre-maxvol [1]
      Procedure :: maxvol252 => mtrx_maxvol252 !Original Dominant-R from [1]
      Procedure :: dominantr => mtrx_dominantr !Faster Dominant-R
      Procedure :: maxvol2r => mtrx_maxvol2r !Greedy column addition from [5]
      Procedure :: maxelement => mtrx_maxelement !Returns position of maximum (in absolute value) element
      Procedure :: fnorm => mtrx_fnorm !Frobenius norm
      Procedure :: cnorm => mtrx_cnorm !Chebyshev norm
      Procedure :: norm2 => mtrx_norm2 !Spectral norm
      Procedure :: trace => mtrx_trace !Trace
      Procedure :: vol => mtrx_vol !Matrix volume
      Procedure :: reshape => mtrx_reshape !Reshapes unfolding matrix of a tensor
      Procedure :: copy => mtrx_copy !Copyies the matrix
      Procedure :: update1v => mtrx_update1v !Rank 1 update with vector type arrays
      Procedure :: unite => mtrx_unite !Unites two matrix arrays into one
      Procedure :: replace => mtrx_replace !Replaces matrix elements according to mask
  End type
  
  !Matrix transpose
  interface operator (.T.)
    module procedure Mtrx_transp
  end interface
  
  !Matrix inverse. Also works with rectangular matrices.
  interface operator (.I.)
    module procedure Mtrx_pinverse
  end interface
  
  !Multiplication by inverse on the right
  interface operator (.dI.)
    module procedure Mtrx_dotinverse
  end interface
  
  !Multiplication by inverse on the left
  interface operator (.Id.)
    module procedure Mtrx_inversedot
  end interface
  
  !Elementwise multiplication
  interface operator (.dot.)
    module procedure Mtrx_dotmul, Mtrx_dotvec, Mtrx_dotvect
  end interface
  
  !Elementwise deletion
  interface operator (.dd.)
    module procedure Mtrx_dotdiv, Mtrx_ddvec, Mtrx_ddvect
  end interface
  
  !Turns vector type or array to matrix type
  interface assignment(=)
    module procedure Mtrx_fromvec, Array_transform
  end interface
  
  !Matrix-matrix, matrix-vector and matrix-number multiplications
  interface operator(*)
    module procedure mtrx_mmtr, mtrx_tmvec, mtrx_mvect, mtrx_tmnum, mtrx_tmnumr, mtrx_mnumt, mtrx_mnumtr
  end interface
  
  !Matrix-number division
  interface operator(/)
    module procedure mtrx_divnum, mtrx_divnumr
  end interface
  
  !Matrix sum
  interface operator(+)
    module procedure mtrx_sum
  end interface
  
  !Matrix subtruction
  interface operator(-)
    module procedure mtrx_sub
  end interface
  
  Contains
  !See brief descriptions above
  
    subroutine mtrx_reshape(this, n, np)
      Class(Mtrx) :: this
      Integer(4), intent(in) :: n
      Integer(4), intent(in), optional :: np
      Integer(4) m
      DOUBLE PRECISION, Allocatable :: d(:,:,:,:)
      if (present(np)) then
        if (n >= np) then
          m = this%n*this%m
          d = reshape(this%d, (/np,n/np,np,m/n/np/), ORDER = (/3,2,1,4/))
          this%d = reshape(d, (/np,m/np/))
          this%n = np
          this%m = m/np
        end if
      else
        m = this%n*this%m/n
        this%d = reshape(this%d, (/n,m/))
        this%n = n
        this%m = m
      end if
    end
  
    !Generates givens rotation
    elemental subroutine givens(c, s, a, b)
      DOUBLE PRECISION, intent(in) :: a, b
      DOUBLE PRECISION, intent(out) :: c, s
      DOUBLE PRECISION tau
      if (b == 0) then
        c = 1
        s = 0
      else if (abs(b) > abs(a)) then
        tau = -a/b
        s = -b/sqrt(b**2 + a**2)
        c = s*tau
      else
        tau = -b/a
        c = a/sqrt(b**2 + a**2)
        s = c*tau
      end if
    end
    
    subroutine Mtrx_fromvec(this, vec)
      Type(Vector), intent(in) :: vec
      Class(Mtrx), intent(out) :: this
      Integer(4) i
      Allocate(this%d(vec%n,1))
      do i = 1, vec%n
        this%d(i, 1) = vec%d(i)
      end do
      this%n = vec%n
      this%m = 1
    end
  
    function mtrx_fnorm(this) Result(res)
      Class(Mtrx) :: this
      DOUBLE PRECISION res, dlange
      DOUBLE PRECISION, dimension(:) :: work(1)
      res = dlange('F', this%n, this%m, this%d, this%n, work)
    end
    
    function mtrx_cnorm(this) Result(res)
      Class(Mtrx) :: this
      DOUBLE PRECISION res, dlange
      DOUBLE PRECISION, dimension(:) :: work(this%n)
      res = dlange('M', this%n, this%m, this%d, this%n, work)
    end
    
    function mtrx_norm2(this) Result(res)
      Class(Mtrx) :: this
      Type(Mtrx) :: x, y, z
      DOUBLE PRECISION res
      call this%svd(x, y, z)
      res = y%d(1,1)
      call x%deinit()
      call y%deinit()
      call z%deinit()
    end
    
    function mtrx_trace(this) Result(res)
      Class(Mtrx) :: this
      DOUBLE PRECISION res
      Integer(4) i
      res = 0
      do i = 1, min(this%n, this%m)
        res = res + this%d(i,i)
      end do
    end
    
    subroutine mtrx_maxelement(this, k, l)
      Class(Mtrx) :: this
      Integer(4) k, l, n , m, i, j
      DOUBLE PRECISION res
      res = abs(this%d(1, 1))
      k = 1
      l = 1
      n = this%n
      m = this%m
      do j = 1, m
        do i = 1, n
          if (abs(this%d(i, j)) > res) then
            res = abs(this%d(i, j))
            k = i
            l = j
          end if
        end do
      end do
    end
    
    subroutine mtrx_cmaxvol(this, per, totsteps, maxsteps, cout, addmat, swapt, ABin)
      Class(Mtrx) :: this
      Type(Mtrx) :: C
      Type(IntVec), optional :: per
      Integer(4) k
      Integer(4), intent(in), optional :: maxsteps
      Integer(4), intent(out), optional :: totsteps
      Type(Mtrx), intent(out), optional :: cout
      Type(Mtrx), optional :: addmat
      Integer(4), intent(in), optional :: swapt
      Type(Mtrx), optional :: ABin
      Type(Vector) CJ, CI
      Integer(4) n, m, i, j, maxperm, curperm
      Integer(4) ij1(2), ij2(2)
      DOUBLE PRECISION CIJ
      n = this%n
      m = this%m
      k = m
      maxperm = k
      if (present(maxsteps)) then
        maxperm = maxsteps
      end if
      if (m > n) then
        print *, "error in cmaxvol"
        print *, m, n
      end if
      if (present(ABin)) then
        if (ABin%n == k) then
          C = .T.ABin
        else
          call C%copy(ABin)
        end if
        C%d(1:k,1:k) = 0.0d0
        do i = 1, k
          C%d(i,i) = 1.0d0
        end do
      else
        C = this .dI. this%subarray(k, k)
      end if
      ij1 = maxloc(C%d)
      ij2 = minloc(C%d)
      if (C%d(ij1(1), ij1(2)) > -C%d(ij2(1),ij2(2))) then
        i = ij1(1)
        j = ij1(2)
      else
        i = ij2(1)
        j = ij2(2)
      end if
      CIJ = C%d(i, j)
      curperm = 0
      do while (abs(CIJ) > 1.0d0)
        if (curperm >= maxperm) then
          exit
        end if
        if (i <= C%m) then
          exit
        end if
        call this%swap(1, i, j)
        if (present(per)) then
          call per%swap(i, j)
        end if
        if (present(addmat)) then
          call addmat%swap(swapt, i, j)
        end if
        CJ = tovec(C%subarray(C%n,j,1,j))
        CI = tovec(C%subarray(i,C%m,i,1))
        CI%d(j) = CI%d(j) - 1.0d0
        call C%update1v(-1.0d0/CIJ, CJ, CI)
        call C%swap(1, i, j)
        ij1 = maxloc(C%d)
        ij2 = minloc(C%d)
        if (C%d(ij1(1), ij1(2)) > -C%d(ij2(1),ij2(2))) then
          i = ij1(1)
          j = ij1(2)
        else
          i = ij2(1)
          j = ij2(2)
        end if
        CIJ = C%d(i, j)
        curperm = curperm + 1
      end do
      if (present(totsteps)) then
        totsteps = curperm
      end if
      if (present(cout)) then
        cout = C
      end if
    end
    
    subroutine mtrx_cmaxvolro(this, ro, per, totsteps, maxsteps, cout, addmat, swapt, ABin)
      Class(Mtrx) :: this
      Double precision, intent(in) :: ro
      Type(Mtrx) :: C
      Type(IntVec), optional :: per
      Integer(4) k
      Integer(4), intent(in), optional :: maxsteps
      Integer(4), intent(out), optional :: totsteps
      Type(Mtrx), intent(out), optional :: cout
      Type(Mtrx), optional :: addmat
      Integer(4), intent(in), optional :: swapt
      Type(Mtrx), optional :: ABin
      Type(Vector) CJ, CI
      Integer(4) n, m, i, j, maxperm, curperm
      Integer(4) ij1(2), ij2(2)
      DOUBLE PRECISION CIJ
      n = this%n
      m = this%m
      k = m
      maxperm = k
      if (present(maxsteps)) then
        maxperm = maxsteps
      end if
      if (m > n) then
        print *, "error in cmaxvolro"
        print *, m, n
      end if
      if (present(ABin)) then
        if (ABin%n == k) then
          C = .T.ABin
        else
          call C%copy(ABin)
        end if
        C%d(1:k,1:k) = 0.0d0
        do i = 1, k
          C%d(i,i) = 1.0d0
        end do
      else
        C = this .dI. this%subarray(k, k)
      end if
      ij1 = maxloc(C%d)
      ij2 = minloc(C%d)
      if (C%d(ij1(1), ij1(2)) > -C%d(ij2(1),ij2(2))) then
        i = ij1(1)
        j = ij1(2)
      else
        i = ij2(1)
        j = ij2(2)
      end if
      CIJ = C%d(i, j)
      curperm = 0
      do while (abs(CIJ) > ro)
        if (curperm >= maxperm) then
          exit
        end if
        if (i <= C%m) then
          exit
        end if
        call this%swap(1, i, j)
        if (present(per)) then
          call per%swap(i, j)
        end if
        if (present(addmat)) then
          call addmat%swap(swapt, i, j)
        end if
        CJ = tovec(C%subarray(C%n,j,1,j))
        CI = tovec(C%subarray(i,C%m,i,1))
        CI%d(j) = CI%d(j) - 1.0d0
        call C%update1v(-1.0d0/CIJ, CJ, CI)
        call C%swap(1, i, j)
        ij1 = maxloc(C%d)
        ij2 = minloc(C%d)
        if (C%d(ij1(1), ij1(2)) > -C%d(ij2(1),ij2(2))) then
          i = ij1(1)
          j = ij1(2)
        else
          i = ij2(1)
          j = ij2(2)
        end if
        CIJ = C%d(i, j)
        curperm = curperm + 1
      end do
      if (present(totsteps)) then
        totsteps = curperm
      end if
      if (present(cout)) then
        cout = C
      end if
    end
    
    subroutine mtrx_fbmaxvol(this, per1, per2, teps, r, curc)
      Class(Mtrx) :: this
      Type(IntVec) :: per1, per2
      Double precision, intent(in) :: teps
      Integer(4), intent(out) :: r
      Type(Mtrx), optional :: curc
      Type(Mtrx) :: err
      Type(Vector) :: u, v
      Integer(4) n, m, k, i, j, i1, ij1(2), ij2(2)
      Double precision cn
      
      Logical rchosen
    
      rchosen = .false.
      n = this%n
      m = this%m
      k = min(n, m)
      err%n = 0
      err%m = 0
      call err%copy(this)
      cn = err%cnorm()
      r = k
      do i1 = 1, k
        ij1 = maxloc(err%d)
        ij2 = minloc(err%d)
        if (err%d(ij1(1), ij1(2)) > -err%d(ij2(1),ij2(2))) then
          i = ij1(1)
          j = ij1(2)
        else
          i = ij2(1)
          j = ij2(2)
        end if
        if ((teps*cn >= abs(err%d(i,j))) .and. (.not. rchosen)) then
          r = i1-1
          rchosen = .true.
        end if
        if ((rchosen) .and. (abs(err%d(i,j)) < eps*max(n,m))) then
          exit
        end if
        u = tovec(err%subarray(n, j, 1, j))
        v = tovec(err%subarray(i, m, i, 1))
        call err%update1v(-1.0d0/err%d(i,j),u,v)
        call per1%swap(i1, i)
        call per2%swap(i1, j)
        if (present(curc)) then
          call curc%swap(1, i1, i)
          call curc%swap(2, i1, j)
        end if
      end do
    end
    
    subroutine mtrx_update1v(this, alpha, x, y)
      Class(Mtrx) :: this
      DOUBLE PRECISION, intent(in) :: alpha
      Type(Vector), intent(in) :: x, y
      call dger(this%n, this%m, alpha, x%d, 1, y%d, 1, this%d, this%n)
    end
    
    subroutine mtrx_hmaxvol2(this, tin, k, l, per, CIN, addmask)
      Class(Mtrx) :: this
      Type(Mtrx) :: C
      Type(IntVec), optional :: per
      Integer(4), optional :: CIN
      Type(Mtrx), intent(in), optional :: addmask
      Type(Vector) :: LB, CI, CJ
      Integer(4), intent(in) :: k, l
      Integer(4), intent(in) :: tin
      Integer(4) t
      Integer(4) n, i1, cr
      Integer(4) ij(1)
      DOUBLE PRECISION tau, ls, alpha
      if (((l > this%n) .and. (tin == 1)) .or. ((l > this%m) .and. (tin .ne. 1))) then
        if (.not. present(CIN)) then
          print *, "error in hmaxvol2"
        end if
      end if
      if (l < k) then
        print *, "error in hmaxvol2"
      end if
      if (tin .eq. 1) then
        t = 1
        n = this%n
        if ((present(CIN)) .and. (CIN == 1)) then
          C = this%subarray(n, k, k+1, 1)
        else
          C = this%subarray(n, k, k+1, 1) .dI. this%subarray(k, k)
        end if
      else
        t = 2
        n= this%m
        if ((present(CIN)) .and. (CIN == 1)) then
          C = .T.this%subarray(k, n, 1, k+1)
        else
          C = .T.(this%subarray(k, k) .Id. this%subarray(k, n, 1, k+1))
        end if
      end if
      LB = (C .dot. C)*evec(k)
      do cr = 1, l-k
        ij = maxloc(LB%d)
        i1 = ij(1)
        call this%swap(t, i1+k, cr+k)
        if (present(per)) then
          call per%swap(i1+k, cr+k)
        end if
        if (present(addmask)) then
          call addmask%swap(tin, i1+k, cr+k)
        end if
        ls = LB%d(i1)
        LB%d(i1) = LB%d(cr)
        LB%d(cr) = 0.0d0
        CI = tovec(C%subarray(i1, k, i1, 2))
        alpha = C%d(i1,1)
        C%d(i1, :) = C%d(cr, :)
        C%d(cr,:) = 0
        call dlarfg(k, alpha, CI%d, 1, tau)
        call CJ%init(k)
        CJ%d(1) = 1.0d0
        CJ%d(2:k) = CI%d
        call C%update1v(-tau, C*CJ, CJ)
        call CJ%deinit()
        CI = tovec(C%subarray(C%n, 1))
        call LB%update1(-ls/(1.0d0+ls),CI .dot. CI)
        ls = 1.0d0/sqrt(1.0d0 + ls)
        C%d(:,1) = C%d(:,1)*ls
      end do
    end
    
    subroutine mtrx_maxvol2(this, t, k, l, per, addmask)
      Class(Mtrx) :: this
      Type(Mtrx) :: A, C, CI, CB, CBI, CS
      Type(IntVec), optional :: per
      Type(Mtrx), intent(in), optional :: addmask
      Type(Vector) :: LB
      Integer(4), intent(in) :: k, l
      Integer(4), intent(in) :: t
      Integer(4) n, m, i, j, i1, j1, cr, tm
      Integer(4) ij(1)
      DOUBLE PRECISION tmp, ls
      Integer(4) tmpi
      n = this%n
      m = this%m
      if (((l > n) .and. (t == 1)) .or. ((l > m) .and. (t .ne. 1))) then
        print *, "error in maxvol2"
      end if
      if (l < k) then
        print *, "error in maxvol2"
      end if
      call LB%init(max(m,n))
      A = this%subarray(k, k)
      if (t .eq. 1) then
      !Добавляем строки
        tm = 1
        C = this%subarray(n, k)
        C = C .dI. A
      else
      !Добавляем столбцы
        n = m
        tm = 2
        C = this%subarray(k, n)
        C = A .Id. C
        C = .T.C
      end if
      do i1 = k+1, n
        LB%d(i1) = 1.0d0
      end do
      do i1 = k+1, n
        do j1 = 1, k
          LB%d(i1) = LB%d(i1) + C%d(i1,j1) ** 2
        end do
      end do
      call CB%init(n, l)
      call dlacpy('A', n, k, C%d, n, CB%d, n)
      do cr = k+1, l
        ij = maxloc(LB%d)
        i = ij(1)
        call this%swap(tm, i, cr)
        if (present(per)) then
          if (tm == 1) then
            tmpi = per%d(i)
            per%d(i) = per%d(cr)
            per%d(cr) = tmpi
          else
            tmpi = per%d(i)
            per%d(i) = per%d(cr)
            per%d(cr) = tmpi
          end if
        end if
        if (present(addmask)) then
          call addmask%swap(tm, i, cr)
        end if
        ls = LB%d(i)
        CI = CB%subarray(i, l, i, 1)
        CBI = (.T.CI)/ls
        CS = CB * CBI
        !CB = CB - (CS * CI)
        call CB%update1v(-1.0d0, tovec(CS), tovec(CI))
        do j = 1, n
          LB%d(j) = LB%d(j) - ls*(CS%d(j,1))**2
        end do
        
        do i1 = 1, n
          CB%d(i1, cr) = CS%d(i1, 1)
        end do
        call CB%swap(1, i, cr)
        tmp = LB%d(cr)
        LB%d(cr) = LB%d(i)
        LB%d(i) = tmp
        LB%d(cr) = 0.0d0
      end do
    end
    
!Old version of dominantc.
!     subroutine mtrx_dominantc(this, t, k, l, nai, perm1, perm2, steps, maxstepsin)
!       Class(Mtrx) :: this
!       Type(Mtrx) :: A, C, CI, CJ, CJ1, CJ2, CI1, CI2, CB, CBI, CBJ, CBJ2
!       Type(Vector), optional :: perm1, perm2
!       Type(Vector) :: LB, LC
!       Integer(4), intent(in), optional :: nai
!       Integer(4), intent(in) :: k, l
!       Integer(4), intent(in) :: t
!       Integer(4), intent(out), optional :: steps
!       Integer(4), intent(in), optional :: maxstepsin
!       Integer(4) n, m, i, j, i1, j1, tm, maxsteps, cursteps, na
!       DOUBLE PRECISION tmp, k11, k12, k21, k22, li, ro, ro1
!       na = 0
!       if (present(nai)) then
!         na = nai
!       end if
!       n = this%n
!       m = this%m
!       if (((l > n) .and. (t == 1)) .or. ((l > m) .and. (t .ne. 1))) then
!         print *, "error in maxvol2.5.0"
!       end if
!       if (l < k) then
!         print *, "error in maxvol2.5.0"
!       end if
!       call LB%init(max(m,n))
!       if (t .eq. 1) then
!       !Change rows
!         tm = 1
!         A = this%subarray(l, k)
!         C = this%subarray(n, k)
!         C = C .dI. A
!       else
!       !Change columns
!         n = m
!         tm = 2
!         A = this%subarray(k, l)
!         C = this%subarray(k, n)
!         C = A .Id. C
!         C = .T.C
!       end if
!       do i1 = l+1, n
!         LB%d(i1) = 1.0d0
!       end do
!       do i1 = l+1, n
!         do j1 = 1, l
!           LB%d(i1) = LB%d(i1) + C%d(i1,j1) ** 2
!         end do
!       end do
!       CB%n = 0
!       CB%m = 0
!       call CB%copy(C)
!       call C%deinit()
!       call LC%init(l)
!       do i1 = 1, l
!         LC%d(i1) = 1.0d0 - CB%d(i1,i1)
!       end do
! 
!       call CB%maxelement(i,j)
!       call LB%maxelement(i1)
!       call LC%maxelement(j1)
!       ro = CB%d(i,j)**2 + LB%d(i)*LC%d(j)
!       if (i <= l) then
!         ro = 0
!       end if
!       ro1 = CB%d(i1,j1)**2 + LB%d(i1)*LC%d(j1)
!       if (ro1 > ro) then
!         ro = ro1
!         i = i1
!         j = j1
!       end if
! 
!       maxsteps = 2*l
!       if (present(maxstepsin)) then
!         maxsteps = maxstepsin
!       end if
!       cursteps = 0
!       do while (ro > 1.0d0)
!         if (cursteps >= maxsteps) exit
!         call this%swap(tm, i, j)
!         if ((present(perm1)) .and. (present(perm2))) then
!           if (tm == 1) then
!             tmp = perm1%d(i)
!             perm1%d(i) = perm1%d(j)
!             perm1%d(j) = tmp
!           else
!             tmp = perm2%d(i)
!             perm2%d(i) = perm2%d(j)
!             perm2%d(j) = tmp
!           end if
!         end if
!         k11 = CB%d(i,j)/ro
!         k12 = (1.0d0 - CB%d(j,j))/ro
!         k21 = LB%d(i)/ro
!         k22 = k11
!         CJ = CB%subarray(j, l, j, 1)
!         CJ%d(1, j) = CJ%d(1,j)-1.0d0
!         CJ1 = CJ*k21
!         CJ2 = CJ*k22
!         CI = CB%subarray(i, l, i, 1)
!         CBI = CB*(.T.CI)
!         CI%d(1, j) = CI%d(1,j)-1.0d0
!         CI1 = CI*k11
!         CI2 = CI*k12
!         CBJ = CB%subarray(n, j, 1, j)
!         li = LB%d(i)
!         CBJ2 = CBJ - CBI * CB%d(i,j) / li
!         do i1 = l+1, n
!           LB%d(i1) = LB%d(i1) - CBI%d(i1,1)**2 / li + CBJ2%d(i1,1)**2 * k21
!         end do
!         call LB%swap(i,j)
!         !CB = CB - CBJ*(CI1-CJ1) - CBI*(CI2+CJ2)
!         call CB%update1(1.0d0,CBJ,CJ1-CI1)
!         call CB%update1(-1.0d0,CBI,CJ2+CI2)
!         call CB%swap(1, i, j)
! 
!         do i1 = 1, l
!           LC%d(i1) = 1.0d0 - CB%d(i1,i1)
!         end do
!         call CB%maxelement(i,j)
!         call LB%maxelement(i1)
!         call LC%maxelement(j1)
!         ro = CB%d(i,j)**2 + LB%d(i)*LC%d(j)
!         if (i <= l) then
!           ro = 0
!         end if
!         ro1 = CB%d(i1,j1)**2 + LB%d(i1)*LC%d(j1)
!         if (ro1 > ro) then
!           ro = ro1
!           i = i1
!           j = j1
!         end if
!         
!         cursteps = cursteps + 1
!       end do
!       if (present(steps)) then
!         steps = cursteps
!       end if
!     end

    subroutine mtrx_dominantc(this, t, k, l, nai, per1, per2, steps, maxstepsin)
      Class(Mtrx) :: this
      Type(Mtrx) :: A, CI, CI1, CB, CBI, B
      Type(IntVec), optional :: per1, per2
      Type(Vector) :: LB, LC
      Integer(4), intent(in), optional :: nai
      Integer(4), intent(in) :: k, l
      Integer(4), intent(in) :: t
      Integer(4), intent(out), optional :: steps
      Integer(4), intent(in), optional :: maxstepsin
      Integer(4) n, m, i, j, i1, j1, tm, maxsteps, cursteps, na
      Integer(4) ij(2)
      DOUBLE PRECISION k11, k12
      Integer(4) tmpi
      na = 0
      if (present(nai)) then
        na = nai
      end if
      n = this%n
      m = this%m
      if (((l > n) .and. (t == 1)) .or. ((l > m) .and. (t .ne. 1))) then
        print *, "error in maxvol2.5.1"
      end if
      if (l < k) then
        print *, "error in maxvol2.5.1"
      end if
      call LB%init(max(m,n))
      if (t .eq. 1) then
      !Work in columns
        tm = 1
        A = this%subarray(l, k)
        CB = this%subarray(n, k)
        CB = CB .dI. A
      else
      !Work in rows
        n = m
        tm = 2
        A = this%subarray(k, l)
        CB = this%subarray(k, n)
        CB = A .Id. CB
        CB = .T.CB
      end if
      do i1 = l+1, n
        LB%d(i1) = 1.0d0
      end do
      do i1 = 1, n
        do j1 = 1, l
          LB%d(i1) = LB%d(i1) + CB%d(i1,j1) ** 2
        end do
      end do
      call LC%init(l)
      do i1 = 1, l
        LC%d(i1) = 1.0d0 - CB%d(i1,i1)
      end do
      B = CB .dot. CB
      call B%update1v(1.0d0,LB,LC)
      ij = maxloc(B%d)
      i = ij(1)
      j = ij(2)
      maxsteps = 2*l
      if (present(maxstepsin)) then
        maxsteps = maxstepsin
      end if
      cursteps = 0
      do while (B%d(i,j) > 1.0d0)
        if (cursteps >= maxsteps) exit
        call this%swap(tm, i, j)
        if ((present(per1)) .and. (present(per2))) then
          if (tm == 1) then
            tmpi = per1%d(i)
            per1%d(i) = per1%d(j)
            per1%d(j) = tmpi
          else
            tmpi = per2%d(i)
            per2%d(i) = per2%d(j)
            per2%d(j) = tmpi
          end if
        end if
        
        k11 = LB%d(i)
        CBI = CB%subarray(i, l, i, 1)/k11
        CI = CB*(.T.CBI)
        do i1 = 1, n
          LB%d(i1) = LB%d(i1) - k11*CI%d(i1,1)**2
        end do
        !CB = CB - CI*CB%subarray(i, l, i, 1)
        call CB%update1v(-1.0d0,tovec(CI),tovec(CB%subarray(i, l, i, 1)))
        call CB%swap(1,i,j)
        call CI%swap(1,i,j)
        call LB%swap(i,j)
        LB%d(i) = LB%d(i)+1.0d0
        LB%d(j) = LB%d(j)-1.0d0
        CI1 = 1.0d0*CI
        do i1 = 1, n
          CI%d(i1,1) = CB%d(i1,j)
          CB%d(i1,j) = CI1%d(i1,1)
        end do
        !k12 = 1.0d0/(1 - (LB%d(i)-1.0d0))
        k12 = 1.0d0/(2 - LB%d(i))
        do i1 = 1, n
          LB%d(i1) = LB%d(i1) + k12*CI%d(i1,1)**2
        end do
        !CB = CB + k12*CI*CB%subarray(i, l, i, 1)
        call CB%update1v(k12,tovec(CI),tovec(CB%subarray(i, l, i, 1)))
        !CB = CB + k12*CI*CB%subarray(j, l, j, 1)
        !call CB%swap(1,i,j)
        
        do i1 = 1, l
          LC%d(i1) = 1.0d0 - CB%d(i1,i1)
        end do
        B = CB .dot. CB
        call B%update1v(1.0d0,LB,LC)
        ij = maxloc(B%d)
        i = ij(1)
        j = ij(2)
        cursteps = cursteps + 1
      end do
      if (present(steps)) then
        steps = cursteps
      end if
    end
    
!Another old version of dominantc.
    !subroutine mtrx_maxvol251(this, t, k, l, nai, perm1, perm2, steps, maxstepsin)
    !  Class(Mtrx) :: this
    !  Type(Mtrx) :: A, CB, B
    !  Type(Vector) :: CI, CJ, CJ1, CJ2, CI1, CI2, CBI, CBJ
    !  Type(Vector), optional :: perm1, perm2
    !  Type(Vector) :: LB, LC
    !  Integer(4), intent(in), optional :: nai
    !  Integer(4), intent(in) :: k, l
    !  Integer(4), intent(in) :: t
    !  Integer(4), intent(out), optional :: steps
    !  Integer(4), intent(in), optional :: maxstepsin
    !  Integer(4) n, m, i, j, i1, j1, tm, maxsteps, cursteps, na
    !  Integer(4) ij(2)
    !  DOUBLE PRECISION tmp, k11, k12, k21, k22, li
    !  na = 0
    !  if (present(nai)) then
    !    na = nai
    !  end if
    !  n = this%n
    !  m = this%m
    !  if (((l > n) .and. (t == 1)) .or. ((l > m) .and. (t .ne. 1))) then
    !    print *, "error in maxvol2.5.1"
    !  end if
    !  if (l < k) then
    !    print *, "error in maxvol2.5.1"
    !  end if
    !  call LB%init(max(m,n))
    !  if (t .eq. 1) then
    !  !Change rows
    !    tm = 1
    !    A = this%subarray(l, k)
    !    CB = this%subarray(n, k)
    !    CB = CB .dI. A
    !  else
    !  !Change columns
    !    n = m
    !    tm = 2
    !    A = this%subarray(k, l)
    !    CB = this%subarray(k, n)
    !    CB = A .Id. CB
    !    CB = .T.CB
    !  end if
    !  LB = (CB .dot. CB) * evec(l)
    !  do i1 = l+1, n
    !    LB%d(i1) = LB%d(i1) + 1.0d0
    !  end do
    !  call LC%init(l)
    !  do i1 = 1, l
    !    LC%d(i1) = 1.0d0 - CB%d(i1,i1)
    !  end do
    !  B = CB .dot. CB
    !  call B%update1v(1.0d0,LB,LC)
    !  ij = maxloc(B%d)
    !  i = ij(1)
    !  j = ij(2)
    !  maxsteps = 2*l*2
    !  if (present(maxstepsin)) then
    !    maxsteps = maxstepsin
    !  end if
    !  cursteps = 0
    !  do while (B%d(i,j) > 1.0d0)
    !    if (cursteps >= maxsteps) exit
    !    call this%swap(tm, i, j)
    !    if ((present(perm1)) .and. (present(perm2))) then
    !      if (tm == 1) then
    !        tmp = perm1%d(i)
    !        perm1%d(i) = perm1%d(j)
    !        perm1%d(j) = tmp
    !      else
    !        tmp = perm2%d(i)
    !        perm2%d(i) = perm2%d(j)
    !        perm2%d(j) = tmp
    !      end if
    !    end if
    !    k11 = CB%d(i,j)/B%d(i,j)
    !    k12 = (1.0d0 - CB%d(j,j))/B%d(i,j)
    !    k21 = LB%d(i)/B%d(i,j)
    !    k22 = k11
    !    CJ = tovec(CB%subarray(j, l, j, 1))
    !    CJ%d(j) = CJ%d(j)-1.0d0
    !    CJ1 = CJ*k21
    !    CJ2 = CJ*k22
    !    CI = tovec(CB%subarray(i, l, i, 1))
    !    CBI = CB*CI
    !    CI%d(j) = CI%d(j)-1.0d0
    !    CI1 = CI*k11
    !    CI2 = CI*k12
    !    CBJ = tovec(CB%subarray(n, j, 1, j))
    !    li = LB%d(i)
    !    call CB%update1v(1.0d0,CBJ,CJ1-CI1)
    !    call CB%update1v(-1.0d0,CBI,CJ2+CI2)
    !    call CB%swap(1, i, j)
    !    call CBJ%update1(-k11/k21, CBI)
    !    call LB%update1(-1.0d0/li, CBI .dot. CBI)
    !    call LB%update1(k21, CBJ .dot. CBJ)
    !    call LB%swap(i,j)
    !    do i1 = 1, l
    !      LC%d(i1) = 1.0d0 - CB%d(i1,i1)
    !    end do
    !    B = CB .dot. CB
    !    call B%update1v(1.0d0,LB,LC)
    !    ij = maxloc(B%d)
    !    i = ij(1)
    !    j = ij(2)
    !    cursteps = cursteps + 1
    !  end do
    !  if (present(steps)) then
    !    steps = cursteps
    !  end if
    !end
    
    subroutine mtrx_premaxvol(this, maxr, per, ABout, cout)
      Class(Mtrx) :: this
      Type(IntVec), intent(inout), optional :: per
      Integer(4), intent(in) :: maxr
      Type(Mtrx), intent(out), optional :: ABout !First maxr rows and columns contain inverse of right triangular part of Ahat
      Type(Vector), intent(out), optional :: cout
      Type(IntVec) piv
      Type(Vector) c, dapr, dapr2, bc, q
      Type(Mtrx) AB
      Integer(4) j, rank
      Integer(4) n, m
      Integer(4) ij(1)
      DOUBLE PRECISION tmp, aa
      Integer(4) tmpi
      
      if (maxr > this%m) then
        print *, "error in premaxvol"
      end if
      
      n = this%n
      m = this%m
      
      call c%init(m)
      call AB%init(maxr,m)
      call piv%init(m)
      do j = 1, m
        dapr = tovec(this%subarray(n, j, 1, j))
        c%d(j) = dapr * dapr
        piv%d(j) = j
      end do
      if (present(per)) then
        piv = per
      end if
      rank = 0
      ij = maxloc(c%d)
      j = ij(1)
      
      do while (rank < maxr)
      
        rank = rank + 1
        
        !Put j in place of rank
        tmpi = piv%d(j)
        piv%d(j) = piv%d(rank)
        piv%d(rank) = tmpi
        tmp = c%d(j)
        c%d(j) = c%d(rank)
        c%d(rank) = tmp
        call AB%swap(2, j, rank)
        call this%swap(2, j, rank)
        aa = sqrt(c%d(rank))
        bc = tovec(this%subarray(n, rank, 1, rank))
        if (rank > 1) then
          dapr = tovec(AB%subarray(rank-1, rank, 1, rank))
          q = (bc - this%subarray(n, rank-1)*dapr)/aa
        else
          q = bc/aa
        end if
        dapr2 = q*this%subarray(n,m,1,1)
        dapr2%d(1:rank) = 0
        
        !Recalculate A^{-1}
        if (rank > 1) then
          AB%d(1:rank-1, rank) = -dapr%d(1:rank-1)/aa
        end if
        AB%d(rank, rank) = 1.0d0/aa
        
        !Recalculate c
        c%d(rank) = 0
        c%d(rank+1:m) = c%d(rank+1:m) - dapr2%d(rank+1:m)**2
        
        !Recalculate AB
        dapr2 = dapr2/aa
        if (rank > 1) then
          call dger(rank-1, AB%m, -1.0d0, dapr%d, 1, dapr2%d, 1, AB%d, AB%n)
        end if
        AB%d(rank,rank+1:m) = dapr2%d(rank+1:m)

        ij = maxloc(c%d)
        j = ij(1)
      end do
      if (present(per)) then
        per = piv
      end if
      if (present(ABout)) then
        ABout = AB
      end if
      if (present(cout)) then
        cout = c
      end if
    end
    
    !Old version of dominantr.
    subroutine mtrx_maxvol252(this, t, rank, l, per, steps, maxministepsin)
      Class(Mtrx) :: this
      Type(IntVec), intent(inout), optional :: per
      Integer(4), intent(in) :: t, rank, l
      Integer(4), intent(out), optional :: steps
      Integer(4), intent(in), optional :: maxministepsin
      Type(IntVec) piv
      Type(Vector) c, w
      Type(Mtrx) mat, mat2, da, dapr1, dapr2, dapr3, AB, u1, u2, u, A, bnew, cnew, q, bc
      Integer(4) j, i, i1, j1
      Integer(4) n, m, maxministeps, curministeps
      Integer(4) ij(2)
      DOUBLE PRECISION tmp, beta, ro, alpha, mu, f, dapr0, aa, bnew0, cnew0
      Integer(4) tmpi
      
      if (rank .eq. l) then
        !Not yet implemented:
        !call this%maxvol(t, rank)
        return
      end if
      
      if (t .eq. 1) then
        da = .T.(this%subarray(this%n, l))
      else
        da = this%subarray(l, this%m)
      end if
      
      n = da%n
      m = da%m
      
      f = 1.0d0
      call c%init(m)
      call w%init(min(m,n))
      call AB%init(min(m,n),m)
      call piv%perm(m)
      if (present(per)) then
        piv = per
      end if
      maxministeps = floor(log(1.0d0*m)/log(2.0d0))+2*l
      if (present(maxministepsin)) then
        maxministeps = maxministepsin
      end if
      
      mat = da%subarray(n, rank)
      call mat%qr(u, A, 1)
      q = u%subarray(n,rank+1)
      do j = 1, m-rank
        mat = da%subarray(n, j+rank, 1, j+rank)
        mat = (.T.mat) * mat
        c%d(j+rank) = mat%d(1,1)
      end do
      da = (.T.u)*da
      do j = 1, m-rank
        mat = da%subarray(rank, j+rank, 1, j+rank)
        mat = (.T.mat) * mat
        c%d(j+rank) = c%d(j+rank) - mat%d(1,1)
      end do
      
      mat = .I.A
      do j = 1, rank
        do i = 1, rank
          w%d(j) = w%d(j) + mat%d(j, i)**2
        end do
      end do
      call mat2%init(min(m,n), rank)
      do i = 1, rank
        do j = 1, rank
          mat2%d(j,i) = mat%d(j,i)
        end do
      end do
      AB = mat2*da%subarray(rank,m)
      call u%deinit()
      call mat2%deinit()
      
        call mat%deinit()
        call mat%init(rank, m-rank)
        do i = 1, rank
          do j = rank+1, m
            mat%d(i,j-rank) = AB%d(i,j)**2 + c%d(j)*w%d(i)
          end do
        end do
        ij = maxloc(mat%d)
        i1 = ij(1)
        j1 = ij(2)
        ro = sqrt(mat%d(i1,j1))
        j1 = j1 + rank
        call mat%deinit()

        curministeps = 0
        do while ((ro > f) .and. (curministeps < maxministeps))
        
          curministeps = curministeps + 1
          
          !Swap i1 and rank
          do i = i1+1, rank
            tmp = w%d(i-1)
            w%d(i-1) = w%d(i)
            w%d(i) = tmp
            tmpi = piv%d(i-1)
            piv%d(i-1) = piv%d(i)
            piv%d(i) = tmpi
            call AB%swap(1,i-1,i)
            call A%swap(2, i-1, i)
            call this%swap(t, i-1, i)
          end do
          do i = i1, rank-1
            call givens(alpha, beta, A%d(i,i), A%d(i+1,i))
            call mat%init(2,2)
            mat%d(1,1) = alpha
            mat%d(1,2) = -beta
            mat%d(2,1) = beta
            mat%d(2,2) = alpha
            mat2 = mat * A%subarray(i+1,rank,i,i)
            do j = i, rank
              A%d(i,j) = mat2%d(1,j-i+1)
              A%d(i+1,j) = mat2%d(2,j-i+1)
            end do
            mat2 = q%subarray(n,i+1,1,i) * (.T.mat)
            do j = 1, n
              q%d(j,i) = mat2%d(j,1)
              q%d(j,i+1) = mat2%d(j,2)
            end do
            call mat2%deinit()
            call mat%deinit()
          end do
          
          !Swap: rank, rank+1, j1 => j1, rank, rank+1
          tmpi = piv%d(j1)
          piv%d(j1) = piv%d(rank+1)
          piv%d(rank+1) = tmpi
          call this%swap(t, j1, rank+1)
          tmp = c%d(j1)
          c%d(j1) = c%d(rank+1)
          c%d(rank+1) = tmp
          aa = sqrt(c%d(rank+1))
          if (t .eq. 1) then
            bc = .T.(this%subarray(rank+1, n, rank+1, 1))
          else
            bc = this%subarray(n, rank+1, 1, rank+1)
          end if
          call AB%swap(2, j1, rank+1)
          dapr1 = A%d(rank,rank)*AB%subarray(rank,m,rank,rank+2)
          dapr0 = A%d(rank,rank)
          if (rank > 1) then
            dapr3 = A%subarray(rank-1,rank,1,rank)
          end if
          mat = A*AB%subarray(rank,rank+1,1,rank+1)
          do i = 1, rank
            A%d(i, rank) = mat%d(i,1)
          end do
          mat2 = (bc - q%subarray(n,rank)*mat)/aa
          do j = 1, n
            q%d(j,rank+1) = mat2%d(j,1)
          end do
          if (t .eq. 1) then
            mat = .T.(this%subarray(m, n, rank+1, 1))
          else
            mat = this%subarray(n, m, 1, rank+1)
          end if
          mat2 = q%subarray(n,rank)
          mat2 = ((.T.bc) - ((.T.bc)*mat2)*(.T.mat2))/aa
          dapr2 = mat2 * mat
          call mat%deinit()
          call mat2%deinit()
          tmpi = piv%d(rank)
          piv%d(rank) = piv%d(rank + 1)
          piv%d(rank + 1) = tmpi
          call this%swap(t, rank, rank + 1)
          
          !Nullifying column rank under the diagonal (one number)
          call givens(alpha, beta, A%d(rank,rank), aa)
          call mat%init(2,2)
          mat%d(1,1) = alpha
          mat%d(1,2) = -beta
          mat%d(2,1) = beta
          mat%d(2,2) = alpha
          mat2 = q%subarray(n,rank+1,1,rank) * (.T.mat)
          do j = 1, n
            q%d(j,rank) = mat2%d(j,1)
            q%d(j,rank+1) = mat2%d(j,2)
          end do
          call mat2%deinit
          call u%init(1, m-rank)
          u%d(1,1) = dapr0
          do i = 2, m-rank
            u%d(1,i) = dapr1%d(1,i-1)
          end do
          call u2%init(1, m-rank)
          do i = 2, m-rank
            u2%d(1,i) = dapr2%d(1,i)
          end do
          mat2 = mat%subarray(2,1) * u + mat%subarray(2,2,1,2) * u2
          cnew0 = mat2%d(2,1)
          cnew = mat2%subarray(2,m-rank,2,2)
          bnew0 = mat2%d(1,1)
          bnew = mat2%subarray(1,m-rank,1,2)
          A%d(rank,rank) = sqrt(aa**2 + A%d(rank,rank)**2)
          call mat2%deinit()
          call u%deinit()
          call u2%deinit()
          call mat%deinit()
          
          !Recalculate c
          c%d(rank+1) = cnew0**2
          do i = rank + 2, m
            c%d(i) = c%d(i) + cnew%d(1,i-rank-1)**2 - dapr2%d(1,i-rank)**2
          end do
          
          !Recalculate w
          if (rank > 1) then
            mat = A%subarray(rank-1,rank-1)
            u = mat%rtsolve(dapr3)
          end if
          w%d(rank) = 1.0d0/A%d(rank,rank)**2
          if (rank > 1) then
            u2 = AB%subarray(rank-1,rank+1,1,rank+1) + AB%d(rank,rank+1)*u
          end if
          do i = 1, rank - 1
            w%d(i) = w%d(i) + w%d(rank)*u2%d(i,1)**2 - u%d(i,1)**2/dapr0**2
          end do
          
          !Recalculate AB
          mu = AB%d(rank,rank+1)
          if (rank > 1) then
            u1 = AB%subarray(rank-1,rank+1,1,rank+1)
          end if
          AB%d(rank, rank+1) = bnew0/A%d(rank, rank)
          do i = rank+2, m
            AB%d(rank, i) = bnew%d(1,i-rank-1)/A%d(rank, rank)
          end do
          if (rank > 1) then
            u2 = (1.0d0 - mu*AB%d(rank,rank+1))*u - AB%d(rank,rank+1)*u1
          end if
          do i = 1, rank-1
            AB%d(i,rank+1) = u2%d(i,1)
          end do
          bnew = bnew/A%d(rank,rank)
          if ((rank > 1) .and. (rank+1 < m)) then
            mat = u*(dapr1/dapr0 - mu*bnew) - u1*bnew
          end if
          do i = rank+2, m
            do j = 1, rank-1
              AB%d(j, i) = AB%d(j,i) + mat%d(j, i-rank-1)
            end do
          end do
          if (rank > 1) then
            call u%deinit()
            call u2%deinit()
            call mat%deinit()
          end if
          
          call mat%init(rank, m-rank)
          do j = rank+1, m
            do i = 1, rank
              mat%d(i,j-rank) = AB%d(i,j)**2 + c%d(j)*w%d(i)
            end do
          end do
          ij = maxloc(mat%d)
          i1 = ij(1)
          j1 = ij(2)
          ro = sqrt(mat%d(i1,j1))
          j1 = j1 + rank
          call mat%deinit()
          
        end do
      if (present(per)) then
        per = piv
      end if
      if (present(steps)) then
        steps = curministeps
      end if
    end
    
    subroutine mtrx_dominantr(this, t, r, l, per, steps, maxstepsin, AIout, ABout, ABin, cin, roin)
      Class(Mtrx) :: this
      Type(IntVec), intent(inout), optional :: per
      Integer(4), intent(in) :: t, r, l
      Integer(4), intent(out), optional :: steps
      Integer(4), intent(in), optional :: maxstepsin
      Type(Mtrx), intent(out), optional :: AIout, ABout
      Type(Mtrx), optional :: ABin
      Type(Vector), optional :: cin
      Double precision, intent(in), optional :: roin
      Type(Mtrx) da, db
      Type(Vector) c, w
      Type(Mtrx) AB, AI, B
      Type(Mtrx) mat
      Type(Mtrx) q
      Integer(4) info
      Integer(4) ij(2), i1, j1, i
      DOUBLE PRECISION ro, maxro, l1, abij, wi
      Integer(4) cursteps, maxsteps
      Type(Vector) c3, c4, c1l1, abj, aig, aii, abi
      Type(Vector) v, vaig
      Double precision anorm
      
      anorm = this%fnorm()
      if (t .eq. 1) then
        da = .T.(this%subarray(r, l))
        db = .T.(this%subarray(this%n, l, r+1, 1))
      else
        da = this%subarray(l, r)
        db = this%subarray(l, this%m, 1, r+1)
      end if
      if (present(roin)) then
        maxro = roin
      else
        maxro = 1.0d0
      end if
      if (present(ABin)) then
        AI = ABin%subarray(r, r)
        AB = ABin%subarray(r, this%m, 1, r+1)
      else
        call da%qr(q, AI)
        B = (.T.q) * db
        AB = AI%rtsolve(B)
        call dtrtri('U', 'N', r, AI%d, r, info)
      end if
      if (present(cin)) then
        call c%init(this%m-r)
        c%d(1:this%m-r) = cin%d(r+1:this%m)
      else
        c = (evec(l) * (db .dot. db)) - (evec(r) * (B .dot. B))
      end if
      c%d(:) = max(0.0d0,c%d(:)-eps*anorm**2)
      
      w = (AI .dot. AI) * evec(r)

      B = (AB .dot. AB)
      call B%update1v(1.0d0,w,c)
      ij = maxloc(B%d)
      i1 = ij(1)
      j1 = ij(2)
      ro = sign(sqrt(B%d(i1,j1)),AB%d(i1,j1))
      cursteps = 0
      maxsteps = 2*r
      if (present(maxstepsin)) then
        maxsteps = maxstepsin
      end if
      
      do while ((abs(ro) > maxro) .and. (cursteps < maxsteps))
        !Initialization of variables
        abij = AB%d(i1,j1)
        l1 = sqrt(c%d(j1))
        aii = tovec(AI%subarray(i1, r, i1, 1))
        wi = w%d(i1)
        aig = AI*(aii/wi)
        abi = tovec(AB%subarray(i1, AB%m, i1, 1))
        abj = tovec(AB%subarray(r, j1, 1, j1))
        c1l1 = (tovec(db%subarray(l, j1, 1, j1)) - da*abj) * db
        c1l1%d(:) = min(c1l1%d(:), l1*sqrt(c%d(:)))
        c1l1%d(:) = max(c1l1%d(:), -l1*sqrt(c%d(:)))
        !Swap columns
        do i = 1, r
          AB%d(i,j1) = 0.0d0
        end do
        AB%d(i1,j1) = 1.0d0
        abi%d(j1) = 1.0d0
        mat = da%subarray(l, i1, 1, i1)
        do i = 1, l
          da%d(i, i1) = db%d(i, j1)
          db%d(i, j1) = mat%d(i, 1)
        end do
        call this%swap(t,i1,j1+r)
        if (present(per)) then
          call per%swap(i1,j1+r)
        end if
        !Recalculate c
        call c3%copy(c1l1)
        c3 = c3 * (1 + abij/ro)
        call c3%update1(-l1**2/ro, abi)
        c3%d(j1) = (1-1/ro)*l1**2
        call c4%copy(c1l1)
        c4 = c4 * (w%d(i1)/(B%d(i1,j1)+abij*ro))
        call c4%update1(1.0d0/ro, abi)
        c4%d(j1) = 1.0d0+1.0d0/ro
        c = c - (c3 .dot. c4)
        c%d(:) = max(0.0d0,c%d(:)-eps*anorm**2)
        
        !Recalculate AI
        call v%copy(abj/ro)
        call v%update1(wi*l1**2/(B%d(i1,j1)+abij*ro), aig)
        v%d(i1) = 1 - 1.0d0/ro
        call AI%update1v(-1.0d0, v, aii)
        !Recalculate w
        w = w + ((wi*v) .dot. (v - 2.0d0 * aig))
        
        !Recalculate AB
        call AB%update1v(-1.0d0, v, abi)
        c3%d(j1) = -l1**2/ro
        call vaig%copy(abj/ro)
        call vaig%update1(-(abij**2 + abij*ro)/(B%d(i1,j1)+abij*ro), aig)
        vaig%d(i1) = 1 - 1.0d0/ro - aig%d(i1)
        call AB%update1v(-wi/(ro+abij), vaig, c3)

        B = (AB .dot. AB)
        call B%update1v(1.0d0,w,c)
        ij = maxloc(B%d)
        i1 = ij(1)
        j1 = ij(2)
        ro = sign(sqrt(B%d(i1,j1)),AB%d(i1,j1))
        cursteps = cursteps + 1
      end do
      if (present(steps)) then
        steps = cursteps
      end if
      
      if (present(AIout)) then
        AIout = AI
      end if
      if (present(ABout)) then
        ABout = AB
      end if
    end
    
    subroutine mtrx_maxvol2r(this, t, k, l)
      Class(Mtrx) :: this
      Type(Mtrx) :: A, C
      Type(Vector) :: LB
      Integer(4), intent(in) :: k, l
      Integer(4), intent(in) :: t
      Integer(4) n, m, i, i1, j1, cr, tm
      n = this%n
      m = this%m
      if ((l > n) .or. (l > m)) then
        print *, "error in maxvol2r"
      end if
      if (l < k) then
        print *, "error in maxvol2r"
      end if
      if (t .eq. 1) then
        tm = 1
      else
        tm = 2
      endif
      do cr = k+1, l
        if (tm .eq. 1) then
        !Add rows
          A = this%subarray(cr-1, l)
          C = this%subarray(n, l)
          C = C .dI. A
        else
        !Add columns
          A = this%subarray(l, cr-1)
          C = this%subarray(l, n)
          C = A .Id. C
          C = .T.C
        end if
        call LB%init(n)
        do i1 = 1, n
          do j1 = 1, k
            LB%d(i1) = LB%d(i1) + C%d(i1,j1) ** 2
          end do
        end do
        call LB%maxelement(i)
        call LB%deinit()
        call this%swap(tm, i, cr)
      end do
    end
    
!Old version of maxvolproj.
    !subroutine mtrx_maxvolproj(this, r, k, l, perm1, perm2, bigsteps, totsteps)
    !  Class(Mtrx) :: this
    !  Integer(4), intent(in) :: r, k, l
    !  Type(Vector), intent(out) :: perm1, perm2
    !  Integer(4), optional :: bigsteps, totsteps
    !  Integer(4) :: cursteps, curministeps, steps, prevsteps
    !  Integer(4) cr, maxst, maxperm, curperm
    !  Integer(4) m, n
    !  Type(Mtrx) cols, rows
    !  Type(Mtrx) ai, u, s, v
    !  n = this%n
    !  m = this%m
    !  if ((k .eq. r) .or. (l .eq. r)) then
    !    call this%maxvolrect(k, l, perm1, perm2, bigsteps, totsteps)
    !    return
    !  end if
    !  call perm1%perm(this%n)
    !  call perm2%perm(this%m)
    !  call this%maxvol(r, perm1, perm2, cursteps, curministeps, 4)
    !  call perm1%deinit()
    !  call perm2%deinit()
    !  call perm1%perm(this%n)
    !  call perm2%perm(this%m)
    !  rows = this%subarray(r, this%m)
    !  call rows%maxvol251(2, r, l, 0, perm1, perm2, curperm, l)
    !  call this%permcols(perm2, 1, l)
    !  call perm1%deinit()
    !  call perm2%deinit()
    !  call perm1%perm(this%n)
    !  call perm2%perm(this%m)
    !  cols = this%subarray(this%n, l) * (.I.this%subarray(r, l))
    !  call cols%maxvol251(1, r, k, 0, perm1, perm2, curperm, k)
    !  call this%permrows(perm1, 1, k)
    !  call perm1%deinit()
    !  call perm2%deinit()
    !  call cols%deinit()
    !  call rows%deinit()
    !  call cols%init(this%n, l)
    !  call rows%init(k, this%m)
    !  cursteps = 1
    !  maxst = 4
    !  if ((present(bigsteps)) .and. (bigsteps > 0)) then
    !    maxst = bigsteps
    !  end if
    !  maxperm = max(k,l)
    !  if ((present(totsteps)) .and. (totsteps > 0)) then
    !    maxperm = totsteps
    !  end if
    !  if ((k > n) .or. (l > m)) then
    !    print *, "error in maxvolproj"
    !  end if
    !  steps = 1
    !  prevsteps = 1
    !  do while (steps > 0)
    !    steps = 0
    !    do cr = 1, 2
    !      curperm = 0
    !      if (cr .eq. 2) then
    !        cols = this%subarray(this%n, l)
    !        ai = this%subarray(k,l)
    !        call ai%svd(u,s,v)
    !        v = v%subarray(r, l)
    !        cols = cols*(.T.v)
    !        call perm1%perm(this%n)
    !        call perm2%perm(this%m)
    !        call cols%maxvol251(1,r,k,0,perm1,perm2,curperm,k)
    !        call this%permrows(perm1, 1, k)
    !        call perm1%deinit()
    !        call perm2%deinit()
    !      else
    !        rows = this%subarray(k, this%m)
    !        ai = this%subarray(k,l)
    !        call ai%svd(u,s,v)
    !        u = u%subarray(k, r)
    !        rows = (.T.u)*rows
    !        call perm1%perm(this%n)
    !        call perm2%perm(this%m)
    !        call rows%maxvol251(2,r,l,0,perm1,perm2,curperm,l)
    !        call this%permcols(perm1, 1, l)
    !        call perm1%deinit()
    !        call perm2%deinit()
    !      end if
    !      steps = steps + curperm
    !      prevsteps = prevsteps + curperm
    !      if ((steps == 0) .and. (cr == 2)) then
    !        exit
    !      end if
    !      if ((prevsteps == 0) .and. (cr == 1)) then
    !        steps = 0
    !        exit
    !      end if
    !      if (cursteps >= maxst) then
    !        steps = 0
    !        exit
    !      end if
    !      if (cr == 1) then
    !        prevsteps = 0
    !      end if
    !      cursteps = cursteps + 1
    !    end do
    !    curministeps = curministeps + steps
    !  end do
    !  if (present(bigsteps)) then
    !    bigsteps = cursteps
    !  end if
    !  if (present(totsteps)) then
    !    totsteps = curministeps
    !  end if
    !end
  
    subroutine mtrx_swap(this, t, a, b)
      Class(Mtrx) :: this
      Integer(4), intent(in) :: a, b
      Integer(4) n, m, i
      Integer(4), intent(in) :: t
      DOUBLE PRECISION tmp
      n = this%n
      m = this%m
      if (t .eq. 1) then
        if ((a > n) .or. (b > n)) then
          print *, "error in swap_rows_mtrx"
          return
        end if
        do i = 1, m
          tmp = this%d(a, i)
          this%d(a, i) = this%d(b, i)
          this%d(b, i) = tmp
        end do
      end if
      if (t .eq. 2) then
        if ((a > m) .or. (b > m)) then
          print *, "error in swap_columns_mtrx"
          return
        end if
        do i = 1, n
          tmp = this%d(i, a)
          this%d(i, a) = this%d(i, b)
          this%d(i, b) = tmp
        end do
      end if
    end
    
    subroutine mtrx_permrows(this, per, t, kin)
      Class(Mtrx) :: this
      Type(Mtrx) mat
      Type(IntVec), intent(in) :: per
      Integer(4) n, m, i, j
      Integer(4), intent(in) :: t
      Integer(4), intent(in), optional :: kin
      Integer(4) k
      Real(8) tmpr
      Integer(4) tmpi
      n = this%n
      m = this%m
      k = n
      if (present(kin)) then
        k = kin
      end if
      if (per%n .ne. n) then
        print *, "error in perm_rows_mtrx", per%n, n
        return
      end if
      call mat%init(n,m)
      do i = 1, k
        tmpi = per%d(i)
        if (t == 1) then
          do j = 1, m
            tmpr = mat%d(i,j)
            mat%d(i,j) = this%d(tmpi,j)
            if (tmpi > k) then
              this%d(tmpi,j) = tmpr
            end if
          end do
        else
          do j = 1, m
            mat%d(tmpi,j) = this%d(i,j)
          end do
        end if
      end do
      do j = 1, m
        do i = 1, k
          this%d(i, j) = mat%d(i, j)
        end do
      end do
      call mat%deinit()
    end
    
    subroutine mtrx_permcols(this, per, t, kin)
      Class(Mtrx) :: this
      Type(Mtrx) mat
      Type(IntVec), intent(in) :: per
      Integer(4) n, m, i, j, k
      Integer(4), intent(in) :: t
      Integer(4), intent(in), optional :: kin
      Real(8) tmpr
      Integer(4) tmpi
      n = this%n
      m = this%m
      k = m
      if (present(kin)) then
        k = kin
      end if
      if (per%n .ne. m) then
        print *, "error in perm_cols_mtrx", per%n, m
        return
      end if
      call mat%init(n,m)
      do i = 1, k
        tmpi = per%d(i)
        if (t == 1) then
          do j = 1, n
            tmpr = mat%d(j,i)
            mat%d(j,i) = this%d(j,tmpi)
            if (tmpi > k) then
              this%d(j, tmpi) = tmpr
            end if
          end do
        else
          do j = 1, n
            mat%d(j,tmpi) = this%d(j,i)
          end do
        end if
      end do
      do j = 1, k
        do i = 1, n
          this%d(i, j) = mat%d(i, j)
        end do
      end do
      call mat%deinit()
    end
    
    !Putting print inside this function breaks it (it does not even start!)
    function mtrx_vol(this, rank) Result(res)
      Class(Mtrx), intent(in) :: this
      Integer(4), intent(in), optional :: rank
      Type(Mtrx) :: u, s, v
      DOUBLE PRECISION res
      Integer(4) n, i
      if (.not. present(rank)) then
        if (this%n > this%m) then !replace with halfqr and halflq
          call this%qr(u, s)
        else
          call this%lq(s,u)
        end if
        n = min(this%n, this%m)
      else
        call this%svd(u,s,v)
        n = rank
      end if
      res = 1
      do i = 1, n
        res = res * s%d(i,i)
      end do
      res = abs(res)
    end
  
    function mtrx_tauinverse(this, tau) Result(res)
      Class(Mtrx), intent(in) :: this
      Type(Mtrx) :: u, s, vt, res
      DOUBLE PRECISION tau
      Integer(4) n, i
      call this%svd(u, s, vt)
      n = min(this%n, this%m)
      do i = 1, n
        if (s%d(i, i) > tau*s%d(1, 1)) then
          s%d(i, i) = 1.0d0 / s%d(i, i)
        else
          s%d(i, i) = 0
        end if
      end do
      res = (.T.vt) * (.T.s) * (.T.u)
    end
  
    recursive function mtrx_pinverse(this) Result(res)
      Class(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Type(Mtrx) q
      Integer(4) n, info
      DOUBLE PRECISION, Allocatable :: work(:), ipiv(:)
      if (this%n < this%m) then
        res = .T.(.I.(.T.this))
      else if (this%n .eq. this%m) then
        n = this%n
        call res%init(n, n)
        Allocate(work(n))
        Allocate(ipiv(n))
        res = 1.0d0*this
        call dgetrf(n, n, res%d, n, ipiv, info)
        if (info /= 0) then
          stop 'Inverse error: matrix is singular!'
        end if
        call dgetri(n, res%d, n, ipiv, work, n, info)
        if (info /= 0) then
          stop 'Inversion failed!'
        end if
        Deallocate(work)
        Deallocate(ipiv)
      else
        call this%qr(q, res)
        res = (.I.res)*(.T.q)
      end if
    end
    
    !Would be good to enable permutations as input
    !Or create another subroutine (same name?) for that
    function mtrx_subarray(this, n, m, k, l) Result(res)
      Class(Mtrx), intent(in) :: this
      Integer(4), intent(in) :: n, m
      Integer(4), intent(in), optional :: k, l
      Integer(4) n1, m1
      Type(Mtrx) res
      if ((.not. present(k)) .or. (.not. present(l))) then
        n1 = 1
        m1 = 1
      else
        n1 = k
        m1 = l
      end if
      res%n = n - n1 + 1
      res%m = m - m1 + 1
      Allocate(res%d(res%n,res%m))
      res%d = this%d(n1:n,m1:m)
    end
    
    function mtrx_subcol(this, cols) Result(res)
      Class(Mtrx), intent(in) :: this
      Type(Vector), intent(in) :: cols
      Type(Mtrx) :: res
      Integer(4) i, j
      call res%init(this%n, cols%n)
      if (cols%n <= this%m) then
        do j = 1, res%m
          do i = 1, res%n
            res%d(i, j) = this%d(i, floor(cols%d(j)+0.5d0))
          end do
        end do
      else
        print *, "error subcolumn_mtrx", this%n, this%m, cols%n
      end if
    end
  
    subroutine mtrx_svd(this, u, s, vt)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: u, s, vt
      DOUBLE PRECISION, Allocatable :: ds(:), da(:,:)
      DOUBLE PRECISION, Allocatable :: work(:)
      Integer(4) Lwork, info
      Integer(4) n, m, i, j
      n = this%n
      m = this%m
      Lwork = 5 * (m + n)
      Allocate(work(Lwork))
      Allocate(da(n, m))
      Allocate(ds(n))
      do j = 1, m
        do i = 1, n
          da(i, j) = this%d(i, j)
        end do
      end do
      call u%init(n, n)
      call s%init(n, m)
      call vt%init(m, m)
      call dgesvd('A', 'A', n, m, da, n, ds, u%d, n, vt%d, m, work, Lwork, info)
      Deallocate(work)
      Deallocate(da)
      do i = 1, min(n, m)
        s%d(i, i) = ds(i)
      end do
      Deallocate(ds)
      if (info .ne. 0) then
        print *, "error svd_mtrx", info
        !call backtrace()
        stop
      end if
    end
    
    function mtrx_rtsolve(this, b) Result(res)
      Class(Mtrx), intent(in) :: this
      Type(Mtrx), intent(in) :: b
      Type(Mtrx) res
      res%n = 0
      res%m = 0
      call res%copy(b)
      call dtrsm('L', 'U', 'N', 'N', b%n, b%m, 1.0d0, this%d, this%n, res%d, b%n)
    end
    
    function mtrx_ltsolve(this, b) Result(res)
      Class(Mtrx), intent(in) :: this
      Type(Mtrx), intent(in) :: b
      Type(Mtrx) res
      res%n = 0
      res%m = 0
      call res%copy(b)
      call dtrsm('L', 'L', 'N', 'N', b%n, b%m, 1.0d0, this%d, this%n, res%d, b%n)
    end
    
    subroutine mtrx_halfqr(this, q, tau, r)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: q
      Type(Mtrx), intent(out) :: r
      Type(Vector), intent(out) :: tau
      DOUBLE PRECISION, Allocatable :: da(:,:)
      DOUBLE PRECISION, Allocatable :: work(:)
      Integer(4) Lwork, info
      Integer(4) n, m
      n = this%n
      m = this%m
      Lwork = 5 * (m + n)
      Allocate(work(Lwork))
      Allocate(da(n, m))
      call tau%init(min(m,n))
      call dlacpy('A', n, m, this%d, n, da, n)
      call dgeqrf(n, m, da, n, tau%d, work, Lwork, info)
      if (info .ne. 0) then
        print *, "error qr(r)_mtrx", info
      end if
      call r%init(min(m, n), m)
      call dlacpy('U', min(n,m), m, da, n, r%d, min(n,m))
      q%n = size(da, 1)
      q%m = size(da, 2)
      q%d = da
      Deallocate(work)
    end
    
    subroutine mtrx_qr(this, q, r, fin)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: q, r
      Integer(4), intent(in), optional :: fin
      Type(Mtrx) qd
      Type(Vector) dtau
      DOUBLE PRECISION, Allocatable :: work(:)
      Integer(4) Lwork, info
      Integer(4) n, m, full
      full = 0
      if (present(fin)) then
        full = 1
      end if
      n = this%n
      m = this%m
      Lwork = 5 * (m + n)
      Allocate(work(Lwork))
      call this%halfqr(qd, dtau, r)
      q%n = n
      q%m = min(m,n)+full*(n-min(m,n))
      Allocate(q%d(q%n,q%m))
      call dlacpy('A', n, min(m,n), qd%d, n, q%d, n)
      call dorgqr(n, min(m,n)+full*(n-min(m,n)), min(m, n), q%d, n, dtau%d, work, Lwork, info)
      Deallocate(work)
      call dtau%deinit()
      if (info .ne. 0) then
        print *, "error qr(q)_mtrx", info
      end if
    end
    
    subroutine mtrx_halflq(this, r, tau, q)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: q, r
      Type(Vector), intent(out) :: tau
      DOUBLE PRECISION, Allocatable :: da(:,:)
      DOUBLE PRECISION, Allocatable :: work(:)
      Integer(4) Lwork, info
      Integer(4) n, m
      Integer(4) i, j
      n = this%n
      m = this%m
      Lwork = 5 * (m + n)
      Allocate(work(Lwork))
      Allocate(da(n, m))
      call tau%init(min(m,n))
      do j = 1, m
        do i = 1, n
          da(i, j) = this%d(i, j)
        end do
      end do
      call dgelqf(n, m, da, n, tau%d, work, Lwork, info)
      if (info .ne. 0) then
        print *, "error lq(l)_mtrx", info
      end if
      call r%init(n, min(m,n))
      do i = 1, n
        do j = 1, min(i,r%m)
          r%d(i,j) = da(i,j)
        end do
      end do
      q = da
      Deallocate(work)
    end
    
    subroutine mtrx_lq(this, r, q)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: q, r
      DOUBLE PRECISION, Allocatable :: work(:)
      Type(Mtrx) qd
      Type(Vector) dtau
      Integer(4) Lwork, info
      Integer(4) n, m
      Integer(4) i, j
      n = this%n
      m = this%m
      Lwork = 5 * (m + n)
      Allocate(work(Lwork))
      call this%halflq(r, dtau, qd)
      call dorglq(min(m,n), m, min(m, n), qd%d, n, dtau%d, work, Lwork, info)
      call q%init(min(m,n), m)
      do j = 1, m
        do i = 1, min(m,n)
          q%d(i, j) = qd%d(i, j)
        end do
      end do
      Deallocate(work)
      call dtau%deinit()
      if (info .ne. 0) then
        print *, "error lq(q)_mtrx", info
      end if
    end
    
    function mtrx_multq(this, q, tau, qp, qf, qtin) Result(res)
      Class(Mtrx), intent(in) :: this
      Type(Mtrx), intent(in) :: q
      Type(Vector), intent(in) :: tau
      Character, intent(in) :: qp, qf
      Character, intent(in), optional :: qtin
      Type(Mtrx) :: res
      DOUBLE PRECISION, Allocatable :: work(:)
      Character qt
      Integer(4) Lwork, info
      if (present(qtin)) then
        qt = 'T'
      else
        qt = 'N'
      endif
      if (qp == 'L') then
        if (qt == 'N') then
          res%n = q%n
        else
          res%n = q%m
        end if
        res%m = this%m
      else
        if (qt == 'N') then
          res%m = q%m
        else
          res%m = q%n
        end if
        res%n = this%n
      end if
      call res%init(res%n, res%m)
      call dlacpy('A', this%n, this%m, this%d, this%n, res%d, res%n)
      Lwork = 5 * (this%m + this%n)
      Allocate(work(Lwork))
      if (qf == 'U') then
        call dormlq(qp, qt, res%n, res%m, tau%n, q%d, q%n, tau%d, res%d, res%n, work, Lwork, info)
      else if (qf == 'D') then
        call dormqr(qp, qt, res%n, res%m, tau%n, q%d, q%n, tau%d, res%d, res%n, work, Lwork, info)
      else
        print *, 'Wrong mult_q parameter (U,D)', qf
        stop
      end if
      Deallocate(work)
      if (info < 0) then
        print *, 'Error in dormlq or dormqr (mult_q)', info
        stop
      end if
    end
    
    subroutine mtrx_hqr(this, q, r, maxr)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: q, r
      Integer(4), intent(in), optional :: maxr
      Type(IntVec) piv
      Type(Vector) c, v, x
      Type(Mtrx) mat, mat2, da
      Integer(4) j, i, rank, info, Lwork
      Integer(4) n, m, maxrank
      DOUBLE PRECISION tmp, beta, tau
      DOUBLE PRECISION, Allocatable :: work(:), dtau(:)
      Integer(4) tmpi
      n = this%n
      m = this%m
      call c%init(m)
      call piv%init(m)
      do j = 1, m
        mat = this%subarray(n, j, 1, j)
        mat = (.T.mat) * mat
        c%d(j) = mat%d(1,1)
        piv%d(j) = j
      end do
      rank = 0
      if (present(maxr)) then
        maxrank = maxr
      else
        maxrank = min(n,m)
      end if
      call c%maxelement(j)
      tau = c%d(j)
      da = 1.0d0*this
      Allocate(dtau(min(m,n)))
      do while ((tau > 0) .and. (rank < maxrank))
      
        rank = rank + 1
        
        tmpi = piv%d(j)
        piv%d(j) = piv%d(rank)
        piv%d(rank) = tmpi
        call da%swap(2, j, rank)
        tmp = c%d(j)
        c%d(j) = c%d(rank)
        c%d(rank) = tmp
        call x%init(n-rank+1)
        mat = da%subarray(n,rank,rank,rank)
        do i = 1,n-rank+1
          x%d(i) = mat%d(i,1)
        end do
        call x%house(v, beta)
        call x%deinit()
        dtau(rank) = beta
        call mat2%init(n-rank+1,1)
        do i = 1, n-rank+1
          mat2%d(i, 1) = v%d(i)
        end do
        call v%deinit()
        mat = da%subarray(n,m,rank,rank)
        !mat = mat - beta*mat2*((.T.mat2)*mat)
        call mat%update1v(-1.0d0*beta,tovec(mat2),tovec((.T.mat2)*mat))
        do j = rank, m
          do i = rank, n
            da%d(i,j) = mat%d(i-rank+1,j-rank+1)
          end do
        end do
        do i = 2, n - rank + 1
          da%d(i+rank-1, rank) = mat2%d(i, 1)
        end do
        call mat2%deinit()
        call mat%deinit()
        do i = rank + 1, m
          c%d(i) = c%d(i) - da%d(rank,i)*da%d(rank,i)
        end do
        
        tau = 0.0d0
        do i = rank+1, m
          if (tau < c%d(i)) then
            j = i
            tau = c%d(i)
          end if
        end do
        
      end do
      call r%init(rank, m)
      do i = 1, rank
        do j = i, m
          r%d(i,j) = da%d(i,j)
        end do
      end do
      call r%permcols(piv, -1)
      Lwork = 5 * (m + n)
      Allocate(work(Lwork))
      call dorgqr(n, rank, rank, da%d, n, dtau, work, Lwork, info)
      call q%init(n, rank)
      do j = 1, rank
        do i = 1, n
          q%d(i, j) = da%d(i, j)
        end do
      end do
      call da%deinit()
      Deallocate(work)
      Deallocate(dtau)
      if (info .ne. 0) then
        print *, "error hqr(q)_mtrx", info
      end if
    end
    
    subroutine mtrx_geqr(this, q, r, maxr, per)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: q, r
      Type(IntVec), intent(inout), optional :: per
      Integer(4), intent(in), optional :: maxr
      Type(IntVec) piv
      Type(Vector) c, w, v, x
      Type(Mtrx) mat, mat2, da, dapr, AB, u1, u2, u
      Integer(4) j, i, rank, i1, j1, i2, j2
      Integer(4) n, m, maxrank, maxministeps, curministeps
      DOUBLE PRECISION tmp, beta, tau, ro, alpha, mu, f
      DOUBLE PRECISION, Allocatable :: dtau(:)
      Integer(4) tmpi
      
      n = this%n
      m = this%m
      
      if (maxr .eq. min(m,n)) then
        call this%qr(q, r)
        return
      end if
      
      if (n .eq. 1) then
        call q%init(1,1)
        q%d(1,1) = 1
        r = 1.0d0 * this
        return
      end if
      f = 2.0d0
      call c%init(m)
      call w%init(min(m,n))
      call AB%init(min(m,n),m)
      call piv%init(m)
      do j = 1, m
        mat = this%subarray(n, j, 1, j)
        mat = (.T.mat) * mat
        c%d(j) = mat%d(1,1)
        piv%d(j) = j
      end do
      if (present(per)) then
        piv = per
      end if
      rank = 0
      if (present(maxr)) then
        maxrank = maxr
      else
        maxrank = min(n,m)
      end if
      maxministeps = floor(log(1.0d0*m)/log(2.0d0))+1
      call c%maxelement(j)
      tau = c%d(j)
      da = 1.0d0*this
      q = eye(n)
      Allocate(dtau(min(m,n)))
      do while ((tau > 0) .and. (rank < maxrank))
      
        rank = rank + 1
        
        tmpi = piv%d(j)
        piv%d(j) = piv%d(rank)
        piv%d(rank) = tmpi
        call da%swap(2, j, rank)
        tmp = c%d(j)
        c%d(j) = c%d(rank)
        c%d(rank) = tmp
        call AB%swap(2, j, rank)
        call x%init(n-rank+1)
        mat = da%subarray(n,rank,rank,rank)
        do i = 1,n-rank+1
          x%d(i) = mat%d(i,1)
        end do
        call x%house(v, beta)
        call x%deinit()
        dtau(rank) = beta
        call mat2%init(n-rank+1,1)
        do i = 1, n-rank+1
          mat2%d(i, 1) = v%d(i)
        end do
        call v%deinit()
        mat = da%subarray(n,m,rank,rank)
        !mat = mat - beta*mat2*((.T.mat2)*mat)
        call mat%update1v(-1.0d0*beta,tovec(mat2),tovec((.T.mat2)*mat))
        do j = rank, m
          do i = rank, n
            da%d(i,j) = mat%d(i-rank+1,j-rank+1)
          end do
        end do
        do i = 2, n - rank + 1
          da%d(i+rank-1, rank) = 0.0d0
        end do
        mat = q%subarray(n,n,1,rank)
        !mat = mat - (mat*(beta*mat2))*(.T.mat2)
        call mat%update1v(-1.0d0*beta,tovec(mat*mat2),tovec(mat2))
        do j = rank, n
          do i = 1, n
            q%d(i,j) = mat%d(i,j-rank+1)
          end do
        end do
        call mat2%deinit()
        call mat%deinit()
        do i = rank + 1, m
          c%d(i) = c%d(i) - da%d(rank,i)*da%d(rank,i)
        end do
        
        w%d(rank) = 1.0d0/da%d(rank,rank)/da%d(rank,rank)
        do i = 1, rank - 1
          w%d(i) = w%d(i) + AB%d(i, rank)*AB%d(i, rank)*w%d(rank)
        end do
        do i = rank+1, m
          AB%d(rank, i) = da%d(rank,i)/da%d(rank, rank)
        end do
        if (rank > 1) then
          mat = AB%subarray(rank-1,rank,1,rank)*da%subarray(rank,m,rank,rank+1)/da%d(rank, rank)
        end if
        do j = 1, rank-1
          do i = rank+1, m
            AB%d(j, i) = AB%d(j,i) - mat%d(j, i-rank)
          end do
        end do
        if (rank > 1) then
          call mat%deinit()
        end if

        mat = AB%subarray(rank,m,1,rank+1)
        call mat%maxelement(i1, j1)
        j1 = j1+rank
        call mat%deinit()
        tau = 0.0d0
        do i = rank+1, m
          if (tau < c%d(i)) then
            j2 = i
            tau = c%d(i)
          end if
        end do
        tau = 0.0d0
        do i = 1, rank
          if (tau < w%d(i)) then
            i2 = i
            tau = w%d(i)
          end if
        end do
        ro = max(c%d(j2)*w%d(i2), abs(AB%d(i1,j1))*abs(AB%d(i1,j1)))
        if (c%d(j2)*w%d(i2) > abs(AB%d(i1,j1))*abs(AB%d(i1,j1))) then
          i1 = i2
          j1 = j2
        end if
        
        curministeps = 0
        do while ((ro > f) .and. (curministeps < maxministeps))
        
          curministeps = curministeps + 1
          j2 = rank+1
          tmpi = piv%d(j1)
          piv%d(j1) = piv%d(j2)
          piv%d(j2) = tmpi
          call da%swap(2, j1, j2)
          tmp = c%d(j1)
          c%d(j1) = c%d(j2)
          c%d(j2) = tmp
          call AB%swap(2, j1, j2)
          do i = i1+1, rank
            tmp = w%d(i-1)
            w%d(i-1) = w%d(i)
            w%d(i) = tmp
            tmpi = piv%d(i-1)
            piv%d(i-1) = piv%d(i)
            piv%d(i) = tmpi
            call AB%swap(1,i-1,i)
            call da%swap(2,i-1,i)
          end do
          do i = i1, rank-1
            call givens(alpha, beta, da%d(i,i), da%d(i+1,i))
            call mat%init(2,2)
            mat%d(1,1) = alpha
            mat%d(1,2) = -beta
            mat%d(2,1) = beta
            mat%d(2,2) = alpha
            mat2 = mat * da%subarray(i+1,m,i,i)
            do j = i, m
              da%d(i,j) = mat2%d(1,j-i+1)
              da%d(i+1,j) = mat2%d(2,j-i+1)
            end do
            mat2 = q%subarray(n,i+1,1,i) * (.T.mat)
            do j = 1, n
              q%d(j,i) = mat2%d(j,1)
              q%d(j,i+1) = mat2%d(j,2)
            end do
            call mat%deinit()
          end do
          if (i1 < rank) then
            call mat2%deinit()
          end if
          i1 = rank
          j1 = rank + 1
          
          rank = rank + 1
          call x%init(n-rank+1)
          mat = da%subarray(n,rank,rank,rank)
          do i = 1,n-rank+1
            x%d(i) = mat%d(i,1)
          end do
          call x%house(v, beta)
          call x%deinit()
          call mat2%init(n-rank+1,1)
          do i = 1, n-rank+1
            mat2%d(i, 1) = v%d(i)
          end do
          call v%deinit()
          mat = da%subarray(n,m,rank,rank)
          !mat = mat - beta*mat2*((.T.mat2)*mat)
          call mat%update1v(-1.0d0*beta,tovec(mat2),tovec((.T.mat2)*mat))
          do j = rank, m
            do i = rank, n
              da%d(i,j) = mat%d(i-rank+1,j-rank+1)
            end do
          end do
          do i = 2, n - rank + 1
            da%d(i+rank-1, rank) = 0.0d0
          end do
          mat = q%subarray(n,n,1,rank)
          !mat = mat - (mat*(beta*mat2))*(.T.mat2)
          call mat%update1v(-1.0d0*beta,tovec(mat*mat2),tovec(mat2))
          do j = rank, n
            do i = 1, n
              q%d(i,j) = mat%d(i,j-rank+1)
            end do
          end do
          call mat2%deinit()
          call mat%deinit()
          rank = rank - 1
          
          tmpi = piv%d(rank)
          piv%d(rank) = piv%d(rank + 1)
          piv%d(rank + 1) = tmpi
          dapr = 1.0d0*da
          call da%swap(2, rank, rank + 1)
          
          call givens(alpha, beta, da%d(rank,rank), da%d(rank+1,rank))
          call mat%init(2,2)
          mat%d(1,1) = alpha
          mat%d(1,2) = -beta
          mat%d(2,1) = beta
          mat%d(2,2) = alpha
          mat2 = mat * da%subarray(rank+1,m,rank,rank)
          do j = rank, m
            da%d(rank,j) = mat2%d(1,j-rank+1)
            da%d(rank+1,j) = mat2%d(2,j-rank+1)
          end do
          mat2 = q%subarray(n,rank+1,1,rank) * (.T.mat)
          do j = 1, n
            q%d(j,rank) = mat2%d(j,1)
            q%d(j,rank+1) = mat2%d(j,2)
          end do
          call mat2%deinit
          call mat%deinit()
          
          c%d(rank+1) = da%d(rank+1,rank+1)**2
          do i = rank + 2, m
            c%d(i) = c%d(i) + da%d(rank+1,i)**2 - dapr%d(rank+1,i)**2
          end do
          
          if (rank > 1) then
            mat = da%subarray(rank-1,rank-1)
            u = mat%rtsolve(dapr%subarray(rank-1,rank,1,rank))
          end if
          w%d(rank) = 1.0d0/da%d(rank,rank)/da%d(rank,rank)
          if (rank > 1) then
            u2 = AB%subarray(rank-1,rank+1,1,rank+1) + AB%d(rank,rank+1)*u
          end if
          do i = 1, rank - 1
            w%d(i) = w%d(i) + u2%d(i,1)**2*w%d(rank) - u%d(i,1)**2/dapr%d(rank,rank)*2
          end do
          
          mu = AB%d(rank,rank+1)
          if (rank > 1) then
            u1 = AB%subarray(rank-1,rank+1,1,rank+1)
          end if
          do i = rank+1, m
            AB%d(rank, i) = da%d(rank,i)/da%d(rank, rank)
          end do
          if (rank > 1) then
            u2 = (1.0d0 - mu*AB%d(rank,rank+1))*u - AB%d(rank,rank+1)*u1
          end if
          do i = 1, rank-1
            AB%d(i,rank+1) = u2%d(i,1)
          end do
          if ((rank > 1) .and. (rank+1 < m)) then
            mat = u*(dapr%subarray(rank,m,rank,rank+2)/dapr%d(rank,rank) - mu*da%subarray(rank,m,rank,rank+2)/da%d(rank,rank))
            !mat = mat - u1*(da%subarray(rank,m,rank,rank+2)/da%d(rank,rank))
            call mat%update1v(-1.0d0/da%d(rank,rank),tovec(u1),tovec(da%subarray(rank,m,rank,rank+2)))
          end if
          do j = 1, rank-1
            do i = rank+2, m
              AB%d(j, i) = AB%d(j,i) + mat%d(j, i-rank-1)
            end do
          end do
          
          mat = AB%subarray(rank,m,1,rank+1)
          call mat%maxelement(i1, j1)
          j1 = j1+rank
          call mat%deinit()
          tau = 0.0d0
          do i = rank+1, m
            if (tau < c%d(i)) then
              j2 = i
              tau = c%d(i)
            end if
          end do
          tau = 0.0d0
          do i = 1, rank
            if (tau < w%d(i)) then
              i2 = i
              tau = w%d(i)
            end if
          end do
          ro = max(c%d(j2)*w%d(i2), abs(AB%d(i1,j1))*abs(AB%d(i1,j1)))
          
          if (c%d(j2)*w%d(i2) > abs(AB%d(i1,j1))*abs(AB%d(i1,j1))) then
            i1 = i2
            j1 = j2
          end if
          
        end do
        
        tau = 0.0d0
        do i = rank+1, m
          if (tau < c%d(i)) then
            j = i
            tau = c%d(i)
          end if
        end do
        
      end do
      call r%init(rank, m)
      do i = 1, rank
        do j = i, m
          r%d(i,j) = da%d(i,j)
        end do
      end do
      call r%permcols(piv, -1)
      q = q%subarray(n, rank)
      call da%deinit()
      Deallocate(dtau)
      if (present(per)) then
        per = piv
      end if
    end
    
    subroutine mtrx_mask(this, n, q)
      Class(Mtrx) :: this
      Integer(4), intent(in) :: n
      DOUBLE PRECISION, intent(in) :: q
      Integer(4) i, j
      DOUBLE PRECISION r1
      call this%init(n, n)
      do j = 1, n
        do i = 1, n
          call random_number(r1)
          if (r1 <= q) then
            this%d(i, j) = 1
          end if
        end do
      end do
    end
    
    subroutine mtrx_badrandom(this, n, m, r)
      Class(Mtrx) :: this
      Integer(4), intent(in) :: n, m
      Integer(4), intent(in), optional :: r
      Integer(4) i, j
      Type(Mtrx) :: U, V
      DOUBLE PRECISION r1
      if (present(r)) then
      call U%init(n, r)
      call V%init(r, m)
      do j = 1, r
        do i = 1, n
          call random_number(r1)
          r1 = r1 * 2 - 1
          U%d(i, j) = r1
        end do
      end do
      do j = 1, m
        do i = 1, r
          call random_number(r1)
          r1 = r1 * 2 - 1
          V%d(i, j) = r1
        enddo
      enddo
      U = U * V
      call V%deinit()
      else
        this%n = n
        this%m = m
        call U%init(n, m)
        do j = 1, m
        do i = 1, n
          call random_number(r1)
          !r1 = r1 * 2 - 1
          U%d(i, j) = r1
        end do
      end do
      end if
      this%n = n
      this%m = m
      this%d = U%d
    end
    
    subroutine mtrx_gauss(this, n, m1)
      Class(Mtrx) :: this
      Type(Vector),Allocatable :: vd(:)
      Integer(4), intent(in) :: n
      Integer(4), intent(in), optional :: m1
      Integer(4) m, i, j
      if (.not. present(m1)) then
        m = n
      else
        m = m1
      endif
      Allocate(vd(m))
      do i = 1, m
        call vd(i)%random(n, 1)
      end do
      call this%init(n, m)
      do j = 1, m
        do i = 1, n
          this%d(i, j) = vd(j)%d(i)
        end do
      end do
    end
  
    subroutine mtrx_random(this, n, m1)
      Class(Mtrx) :: this
      Type(Mtrx) vdd, r
      Type(Vector),Allocatable :: vd(:)
      Integer(4), intent(in) :: n
      Integer(4), intent(in), optional :: m1
      Integer(4) m, i, j
      if (.not. present(m1)) then
        m = n
      else
        m = m1
      endif
      Allocate(vd(m))
      do i = 1, m
        call vd(i)%random(n)
      end do
      call vdd%init(n, m)
      do j = 1, m
        do i = 1, n
          vdd%d(i, j) = vd(j)%d(i)
        end do
      end do
      call vdd%qr(this, r)
      do i = 1, m
        call vd(i)%deinit()
      end do
      Deallocate(vd)
    end
    
    function Mtrx_tovec(this) Result(res)
      Class(Mtrx) :: this
      Type(vector) :: res
      call res%init(this%n*this%m)
      res%d = reshape(this%d, (/ this%n*this%m /))
    end
    
    subroutine Mtrx_copy(this, mat)
      Type(Mtrx), intent(in) :: mat
      Class(Mtrx) :: this
      if (.not. allocated(this%d)) allocate(this%d(mat%n, mat%m))
      if ((this%n .ne. mat%n) .or. (this%m .ne. mat%m)) then
        Deallocate(this%d)
        Allocate(this%d(mat%n, mat%m))
      end if
      this%n = mat%n
      this%m = mat%m
      call dlacpy('A', mat%n, mat%m, mat%d, mat%n, this%d, this%n)
    end
    
    subroutine Array_transform(this, array)
      DOUBLE PRECISION, dimension(:,:), intent(in) :: array
      Class(Mtrx), intent(out) :: this
      this%n = size(array,1)
      this%m = size(array,2)
      this%d = array
    end
  
    subroutine Mtrx_constructor(this, n, m)
      Class(Mtrx) :: this
      Integer(4) n
      Integer(4) m
      this%n = n
      this%m = m
      if (allocated(this%d)) then
        Deallocate(this%d)
      end if
      Allocate(this%d(n, m))
      this%d = 0
    end
    
    subroutine Mtrx_set(this, d)
      Class(Mtrx) :: this
      DOUBLE PRECISION :: d(:, :)
      this%d = d
    end
    
    subroutine Mtrx_destructor(this)
      Class(Mtrx) :: this
      this%n = 0
      this%m = 0
      !if (.not. allocated(this%d)) then
      !  call backtrace()
      !end if
      Deallocate(this%d)
    end
    
    function mtrx_unite(this, t, m2) Result(res)
      Class(Mtrx), intent(in) :: this
      Integer(4), intent(in) :: t
      Type(Mtrx), intent(in) :: m2
      Type(Mtrx) :: res
      Integer(4) i, j
      if ((t == 1) .and. (this%m == m2%m)) then
        Allocate(res%d(this%n+m2%n,this%m))
        call dlacpy('A', this%n, this%m, this%d, this%n, res%d, this%n+m2%n)
        res%n = this%n + m2%n
        res%m = this%m
        do j = 1, m2%m
          do i = 1, m2%n
            res%d(i+this%n,j) = m2%d(i,j)
          end do
        end do
      else if ((t == 2) .and. (this%n == m2%n)) then
        Allocate(res%d(this%n,this%m+m2%m))
        call dlacpy('A', this%n, this%m, this%d, this%n, res%d, this%n)
        res%n = this%n
        res%m = this%m + m2%m
        do j = 1, m2%m
          do i = 1, m2%n
            res%d(i,j+this%m) = m2%d(i,j)
          end do
        end do
      else
        print *, "error mtrx_unite", t, this%n, this%m, m2%n, m2%m
      endif
    end
    
    function mtrx_dotinverse(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res, q, r
      res%n = this%n
      res%m = this%m
      if (this%m == m2%m) then
        if (m2%n > m2%m) then
          call m2%qr(q,r)
          r = .T.r
          res = (.T.(r%ltsolve(.T.this)))*(.T.q)
        else
          call m2%lq(r,q)
          r = .T.r
          res = .T.(r%rtsolve(q*(.T.this)))
        end if
      else
        print *, "error mul_inverse_mtrx"
      endif
    end
    
    function mtrx_inversedot(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res, q, r
      res%n = this%n
      res%m = this%m
      if (this%n == m2%n) then
        if (this%n > this%m) then
          call this%qr(q,r)
          res = r%rtsolve((.T.q)*m2)
        else
          call this%lq(r,q)
          res = (.T.q)*r%ltsolve(m2)
        end if
      else
        print *, "error inverse_mul_mtrx"
      endif
    end
    
    subroutine mtrx_replace(this, m2, coef, t, r, per, per2)
      Class(Mtrx) :: this
      Type(Mtrx) :: m2
      DOUBLE PRECISION, intent(in) :: coef
      Integer(4), intent(in) :: t
      Integer(4), intent(in) :: r
      Type(IntVec), intent(in) :: per, per2
      DOUBLE PRECISION coef2
      Integer(4) i, j, m, tmp, tmp2
      if (t .eq. 1) then
        m = this%m
      else
        m = this%n
      end if
      if (coef .eq. 1.0d0) then
        do i = 1, r
          tmp = per%d(i)
          if (t == 1) then
            do j = 1, m
              tmp2 = per2%d(j)
              if (m2%d(tmp,tmp2) .ne. 0.0d0) then
                this%d(i,j) = m2%d(tmp,tmp2)
              end if
            end do
          else
            do j = 1, m
              tmp2 = per2%d(j)
              if (m2%d(tmp2,tmp) .ne. 0.0d0) then
                this%d(j,i) = m2%d(tmp2,tmp)
              end if
            end do
          end if
        end do
      else
        coef2 = coef - 1.0d0
        do i = 1, r
          tmp = per%d(i)
          if (t == 1) then
            do j = 1, m
              tmp2 = per2%d(j)
              if (m2%d(tmp,tmp2) .ne. 0.0d0) then
                this%d(i,j) = coef * m2%d(tmp,tmp2) - coef2 * this%d(i,j)
              end if
            end do
          else
            do j = 1, m
              tmp2 = per2%d(j)
              if (m2%d(tmp2,tmp) .ne. 0.0d0) then
                this%d(j,i) = coef * m2%d(tmp2,tmp) - coef2 * this%d(j,i)
              end if
            end do
          end if
        end do
      end if
    end
    
    elemental function mtrx_tmnum(this, num) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      DOUBLE PRECISION, intent(in) :: num
      !Integer(4) info
      !call res%copy(this)
      !if (num*0.0d0 .ne. 0.0d0) then
      !  print *, num
      !end if
      !call dlascl('G', 0, 0, 1.0d0, num, this%n, this%m, this%d, this%n, info)
      !if (info .ne. 0) then
      !  print *, "Error in mtrx_tmnum", info
      !end if
      
      !res%d = res%d - (num*this%d)
      !print *, res%fnorm()
      !Wrong value 50/50. And NaN can happen.
      
      res%n = this%n
      res%m = this%m
      res%d = num * this%d
    end
    
    elemental function mtrx_tmnumr(this, numin) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Real(4), intent(in) :: numin
      DOUBLE PRECISION num
      num = numin
      res%n = this%n
      res%m = this%m
      res%d = num * this%d
    end
    
    elemental function mtrx_mnumt(num, this) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      DOUBLE PRECISION, intent(in) :: num
      res%n = this%n
      res%m = this%m
      res%d = num * this%d
    end
    
    elemental function mtrx_mnumtr(num, this) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Real(4), intent(in) :: num
      res%n = this%n
      res%m = this%m
      res%d = num * this%d
    end
    
    elemental function mtrx_divnum(this, num) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      DOUBLE PRECISION, intent(in) :: num
      res = this * (1.0d0 / num)
    end
    
    elemental function mtrx_divnumr(this, num) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Real(4), intent(in) :: num
      res = this * (1.0d0 / num)
    end
    
    elemental function mtrx_dotmul(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res
      res%n = this%n
      res%m = this%m
      !if ((this%m == m2%m) .and. (this%n == m2%n)) then
        res%d = this%d * m2%d
      !else
      !  print *, "error dotmul_mtrx"
      !endif
    end
    
    elemental function mtrx_dotvec(this, v) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Vector), intent(in) :: v
      Type(Mtrx) :: res
      Integer(4) i
      res%n = this%n
      res%m = this%m
      Allocate(res%d(res%n, res%m))
      !if (this%n == v%n) then
        do i = 1, this%m
          res%d(:,i) = this%d(:,i) * v%d(:)
        end do
      !else
      !  print *, "error dotvec_mtrx"
      !endif
    end
    
    !UNTESTED: FAILS ON SVD
    elemental function mtrx_dotvect(v, this) Result(res)
      Type(Vector), intent(in) :: v
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Integer(4) i
      res%n = this%n
      res%m = this%m
      Allocate(res%d(res%n, res%m))
      !if (this%n == v%n) then
        do i = 1, this%n
          res%d(i,:) = this%d(i,:) * v%d(:)
        end do
      !else
      !  print *, "error dotvec_mtrx"
      !endif
    end
    
    elemental function mtrx_ddvec(this, v) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Vector), intent(in) :: v
      Type(Mtrx) :: res
      Integer(4) i
      res%n = this%n
      res%m = this%m
      Allocate(res%d(res%n, res%m))
      !if (this%n == v%n) then
        do i = 1, this%m
          res%d(:,i) = this%d(:,i) / v%d(:)
        end do
      !else
      !  print *, "error ddvec_mtrx"
      !endif
    end
    
    !UNTESTED
    elemental function mtrx_ddvect(v, this) Result(res)
      Type(Vector), intent(in) :: v
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Integer(4) i
      res%n = this%n
      res%m = this%m
      Allocate(res%d(res%n, res%m))
      !if (this%n == v%n) then
        do i = 1, this%n
          res%d(i,:) = this%d(i,:) / v%d(:)
        end do
      !else
      !  print *, "error ddvec_mtrx"
      !endif
    end
    
    elemental function mtrx_dotdiv(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res
      res%n = this%n
      res%m = this%m
      !if ((this%m == m2%m) .and. (this%n == m2%n)) then
        res%d = this%d / m2%d
      !else
      !  print *, "error dotdiv_mtrx"
      !endif
    end
  
    !Add .dT. for transpose and check for size to do matrix-vector multiplication.
    function mtrx_mmtr(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx), intent(in) :: m2
      Type(Mtrx) :: res
      res%n = this%n
      res%m = m2%m
      !if (this%m == m2%n) then
        !res%d = matmul(this%d, m2%d)
        Allocate(res%d(this%n, m2%m))
        call dgemm('N', 'N', this%n, m2%m, this%m, 1.0d0, this%d, this%n, m2%d, m2%n, 0.0d0, res%d, this%n)
      !else
      !  print *, "error mul_mtrx", this%m, m2%n
      !endif
    end
    
    function mtrx_tmvec(this, v) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Vector) :: res
      Type(Vector), intent(in) :: v
      res%n = this%n
      !if (this%m == v%n) then
        !res%d = matmul(this%d, v%d)
        Allocate(res%d(this%n))
        call dgemv('N', this%n, this%m, 1.0d0, this%d, this%n, v%d, 1, 0.0d0, res%d, 1)
      !else
      !  print *, "error mul_mvec"
      !endif
    end
    
    function mtrx_mvect(v, this) Result(res)
      Type(Vector) :: res
      Type(Vector), intent(in) :: v
      Type(Mtrx), intent(in) :: this
      res%n = this%m
      !if (this%n == v%n) then
        Allocate(res%d(this%m))
        call dgemv('T', this%n, this%m, 1.0d0, this%d, this%n, v%d, 1, 0.0d0, res%d, 1)
      !else
      !  print *, "error mul_mvec"
      !endif
    end
    
    elemental function mtrx_sum(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res
      !if ((this%m == m2%m) .and. (this%n == m2%n)) then
        res%n = this%n
        res%m = this%m
        res%d = this%d + m2%d
      !else
      !  print *, "error sum_mtrx"
      !endif
    end
    
    elemental function mtrx_transp(this) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Integer(4) i, j
      res%n = this%m
      res%m = this%n
      Allocate(res%d(res%n, res%m))
      !if (size(this%d,1) .ne. this%n) then
      !  print *, '1', size(this%d,1), this%n
      !end if
      !if (size(this%d,2) .ne. this%m) then
      !  print *, '2', size(this%d,2), this%m
      !end if
      do i = 1, this%n
        do j = 1, this%m
          res%d(j,i) = this%d(i,j)
        end do
      end do
      !res%d = transpose(this%d(:,:this%m))
    end
    
    elemental function mtrx_sub(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res
      !if ((this%m == m2%m) .and. (this%n == m2%n)) then
        res%n = this%n
        res%m = this%m
        res%d = this%d - m2%d
      !else
      !  print *, "error sub_mtrx"
        !call backtrace()
      !  stop
      !endif
    end
    
    elemental function tovec(this) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(vector) :: res
      res%n = this%n*this%m
      Allocate(res%d(res%n))
      !call res%init(this%n*this%m)
      res%d = reshape(this%d, (/ this%n*this%m /))
    end
    
    !UNTESTED: FAILS ON SVD
    elemental function diag(this) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(vector) :: res
      Integer(4) i
      res%n = min(this%n,this%m)
      Allocate(res%d(res%n))
      do i = 1, res%n
        res%d(i) = this%d(i,i)
      end do
    end
    
    !Generates identity matrix (or rectangular with 1 on diagonal)
    elemental function eye(n, m1) Result(res)
      Integer(4), intent(in) :: n
      Integer(4), intent(in), optional :: m1
      Type(Mtrx) :: res
      Integer(4) m, i
      if (present(m1)) then
        m = m1
      else
        m = n
      end if
      res%n = n
      res%m = m
      Allocate(res%d(n,m))
      res%d = 0
      do i = 1, min(m,n)
        res%d(i,i) = 1.0d0
      end do
    end
end
