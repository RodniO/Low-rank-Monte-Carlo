Module ModMtrx
  Use ModVec
  
  !Добавить A_r и A_r^+

  Type Mtrx
    Integer(4) n, m
    DOUBLE PRECISION,Allocatable :: d(:,:)
    Contains
      Procedure :: init => Mtrx_constructor
      Procedure :: deinit => Mtrx_destructor
      Procedure :: set => Mtrx_set
      Procedure :: mtovec => Mtrx_tovec
      Procedure :: gauss => mtrx_gauss
      Procedure :: random => mtrx_random
      Procedure :: mask => mtrx_mask
      Procedure :: badrandom => mtrx_badrandom
      Procedure :: submatrix => mtrx_subm
      Procedure :: subarray => mtrx_subarray
      Procedure :: subcol => mtrx_subcol
      Procedure :: svd => mtrx_svd
      Procedure :: rtsolve => mtrx_rtsolve
      Procedure :: ltsolve => mtrx_ltsolve
      Procedure :: halfqr => mtrx_halfqr
      Procedure :: qr => mtrx_qr
      Procedure :: halflq => mtrx_halflq
      Procedure :: lq => mtrx_lq
      Procedure :: multq => mtrx_multq
      Procedure :: hqr => mtrx_hqr
      Procedure :: geqr => mtrx_geqr
      Procedure :: tauinverse => mtrx_tauinverse
      Procedure :: swap => mtrx_swap
      Procedure :: pertrows => mtrx_pertrows
      Procedure :: pertcols => mtrx_pertcols
      Procedure :: cmaxvol => mtrx_cmaxvol
      Procedure :: fbmaxvol => mtrx_fbmaxvol
      Procedure :: maxvol => mtrx_maxvol
      Procedure :: maxvol2 => mtrx_maxvol2
      Procedure :: hmaxvol2 => mtrx_hmaxvol2
      Procedure :: maxvol25 => mtrx_maxvol25
      Procedure :: dominantc => mtrx_dominantc
      Procedure :: maxvol251 => mtrx_maxvol251
      Procedure :: premaxvol => mtrx_premaxvol
      Procedure :: maxvol252 => mtrx_maxvol252
      Procedure :: dominantr => mtrx_dominantr
      Procedure :: maxvol2r => mtrx_maxvol2r
      Procedure :: maxvolrect => mtrx_maxvolrect
      Procedure :: maxvolproj => mtrx_maxvolproj
      Procedure :: maxelement => mtrx_maxelement
      Procedure :: fnorm => mtrx_fnorm
      Procedure :: cnorm => mtrx_cnorm
      Procedure :: norm2 => mtrx_norm2
      Procedure :: vol => mtrx_vol
      Procedure :: reshape => mtrx_reshape
      Procedure :: copy => mtrx_copy
      Procedure :: tovec => Mtrx_transform
      Procedure :: update1 => mtrx_update1
      Procedure :: update1v => mtrx_update1v
      Procedure :: unite => mtrx_unite
      Procedure :: replace => mtrx_replace
  End type
  
  interface operator (.T.)
    module procedure Mtrx_transp
  end interface
  
  interface operator (.I.)
    module procedure Mtrx_pinverse
  end interface
  
  interface operator (.dI.)
    module procedure Mtrx_dotinverse
  end interface
  
  interface operator (.Id.)
    module procedure Mtrx_inversedot
  end interface
  
  interface operator (.dot.)
    module procedure Mtrx_dotmul, Mtrx_dotvec
  end interface
  
  interface operator (.dd.)
    module procedure Mtrx_dotdiv, Mtrx_ddvec
  end interface
  
  interface assignment(=)
    module procedure Mtrx_fromvec, Array_transform
  end interface
  
  interface operator(*)
    module procedure mtrx_mmtr, mtrx_tmvec, mtrx_mvect, mtrx_tmnum, mtrx_tmnumr, mtrx_mnumt, mtrx_mnumtr
  end interface
  
  interface operator(/)
    module procedure mtrx_divnum, mtrx_divnumr
  end interface
  
  interface operator(+)
    module procedure mtrx_sum
  end interface
  
  interface operator(-)
    module procedure mtrx_sub
  end interface
  
  Contains
  
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
  
    subroutine givens(c, s, a, b)
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
    
    subroutine mtrx_cmaxvol(this, pert, totsteps, maxsteps, cout, addmat, swapt)
      Class(Mtrx) :: this
      Type(Mtrx) :: C
      Type(Vector), optional :: pert
      Integer(4) k
      Integer(4), intent(in), optional :: maxsteps
      Integer(4), intent(out), optional :: totsteps
      Type(Mtrx), intent(out), optional :: cout
      Type(Mtrx), optional :: addmat
      Integer(4), intent(in), optional :: swapt
      Type(Vector) CJ, CI
      Integer(4) n, m, i, j, maxpert, curpert
      Integer(4) ij1(2), ij2(2)
      DOUBLE PRECISION CIJ, tmp
      n = this%n
      m = this%m
      k = m
      maxpert = k
      if (present(maxsteps)) then
        maxpert = maxsteps
      end if
      if (m > n) then
        print *, "error in cmaxvol"
      end if
      C = this .dI. this%subarray(k, k)
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
      curpert = 0
      do while (abs(CIJ) > 1.0d0)
        if (curpert >= maxpert) then
          exit
        end if
        if (i <= C%m) then
          exit
        end if
        call this%swap(1, i, j)
        if (present(pert)) then
          call pert%swap(i, j)
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
        curpert = curpert + 1
      end do
      if (present(totsteps)) then
        totsteps = curpert
      end if
      if (present(cout)) then
        cout = C
      end if
    end
    
    subroutine mtrx_fbmaxvol(this, pert1, pert2, teps, r, curc)
      Class(Mtrx) :: this
      Type(Vector) :: pert1, pert2
      Double precision, intent(in) :: teps
      Integer(4), intent(out) :: r
      Type(Mtrx), optional :: curc
      Type(Mtrx) :: err
      Type(Vector) :: u, v
      Integer(4) n, m, k, i, j, i1, ij1(2), ij2(2)
      Double precision cn
    
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
        if (teps*cn > abs(err%d(i,j))) then
          r = i1-1
          exit
        end if
        u = tovec(err%subarray(n, j, 1, j))
        v = tovec(err%subarray(i, m, i, 1))
        call err%update1v(-1.0d0/err%d(i,j),u,v)
        call pert1%swap(i1, i)
        call pert2%swap(i1, j)
        if (present(curc)) then
          call curc%swap(1, i1, i)
          call curc%swap(2, i1, j)
        end if
      end do
    end
    
    subroutine mtrx_update1(this, alpha, xin, yin)
      Class(Mtrx) :: this
      DOUBLE PRECISION, intent(in) :: alpha
      Type(Mtrx), intent(in) :: xin, yin
      Type(vector) :: x, y
      call xin%tovec(x)
      call yin%tovec(y)
      call dger(this%n, this%m, alpha, x%d, 1, y%d, 1, this%d, this%n)
    end
    
    subroutine mtrx_update1v(this, alpha, x, y)
      Class(Mtrx) :: this
      DOUBLE PRECISION, intent(in) :: alpha
      Type(Vector), intent(in) :: x, y
      call dger(this%n, this%m, alpha, x%d, 1, y%d, 1, this%d, this%n)
    end
  
    subroutine mtrx_maxvol(this, k, pert1, pert2, bigsteps, totsteps, maxsteps, upgraded)
      Class(Mtrx) :: this
      Type(Mtrx) :: A, C, CJ, CI, CC(2)
      Type(Vector), optional :: pert1, pert2
      Type(Vector) :: CC0(2), p(2)
      Integer(4), intent(in) :: k
      Integer(4), intent(in), optional :: maxsteps
      Integer(4), intent(in), optional :: upgraded
      Integer(4), intent(out), optional :: bigsteps, totsteps
      Integer(4) :: cursteps, curministeps, steps, prevsteps
      Integer(4) n, m, i, j, cr, maxst, maxpert, curpert, i1, curi1, maxi1
      Integer(4) dops(2), j1, i20, i21, j20, j21, i2, j2
      Integer(4) ij1(2), ij2(2)
      DOUBLE PRECISION CIJ, tmp, curmaxvol
      n = this%n
      m = this%m
      cursteps = 1
      curministeps = 0
      maxst = 20
      maxpert = k
      if (present(maxsteps)) then
        maxst = maxsteps
      end if
      if ((k > n) .or. (k > m)) then
        print *, "error in maxvol"
      end if
      A = this%submatrix(k, k)
      C = this%submatrix(n, k)
      C = C .dI. A
      do cr = 1, 2
        call CC0(cr)%init(max(m,n))
      end do
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
      steps = 1
      prevsteps = 1
      do while (steps > 0)
        steps = 0
        do cr = 1, 2
          curpert = 0
          !Не работает? Результат хуже: см. unitcrossdet
          !Другой upgraded глючит при ранге 2: в p -1
          if (present(upgraded)) then
            curmaxvol = 1
            curi1 = 1
            if (cr == 1) then
              maxi1 = n/k
            else
              maxi1 = m/k
            end if
            do i1 = 2, maxi1
              A = C%submatrix(i1*k,k,(i1-1)*k+1,1)
              if (A%vol() > curmaxvol) then
                curi1 = i1
                curmaxvol = A%vol()
              end if
            end do
            if (curi1 > 1) then
              do i1 = 1, k
                call this%swap(cr, i1, i1+(curi1-1)*k)
                if (cr == 1) then
                  tmp = pert1%d(i1)
                  pert1%d(i1) = pert1%d(i1+(curi1-1)*k)
                  pert1%d(i1+(curi1-1)*k) = tmp
                else
                  tmp = pert2%d(i1)
                  pert2%d(i1) = pert2%d(i1+(curi1-1)*k)
                  pert2%d(i1+(curi1-1)*k) = tmp
                end if
              end do
              A = this%submatrix(k, k)
              C = this%submatrix(n, k)
              C = C .dI. A
              call C%maxelement(i, j)
              CIJ = C%d(i, j)
            end if
          end if
          do while (abs(CIJ) > 1.0d0)
            if (i == j) then
              exit
            end if
            call this%swap(cr, i, j)
            if ((present(pert1)) .and. (present(pert2))) then
              if (cr == 1) then
                tmp = pert1%d(i)
                pert1%d(i) = pert1%d(j)
                pert1%d(j) = tmp
              else
                tmp = pert2%d(i)
                pert2%d(i) = pert2%d(j)
                pert2%d(j) = tmp
              end if
            end if
            CJ = C%submatrix(C%n,j,1,j)
            CI = C%submatrix(i,C%m,i,1)
            CI%d(1, j) = CI%d(1, j) - 1.0d0
            call C%update1(-1.0d0/CIJ, CJ, CI)
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
            curpert = curpert + 1
            if (curpert >= maxpert) then
              exit
            end if
          end do
          steps = steps + curpert
          prevsteps = prevsteps + curpert
          CC(cr) = C
          do i1 = k+1, C%n
            CI = C%submatrix(i1,k,i1,1)
            CC0(cr)%d(i1) = CI%cnorm()
          end do
          if ((steps == 0) .and. (cr == 2)) then
            exit
          end if
          if ((prevsteps == 0) .and. (cr == 1)) then
            exit
          end if
          if (cursteps >= maxst) then
            exit
          end if
          if (cr == 1) then
            prevsteps = 0
          end if
          A = this%submatrix(k, k)
          if (cr .eq. 2) then
            C = this%submatrix(n, k)
            C = C .dI. A
          else
            C = this%submatrix(k, m)
            C = A .Id. C
            C = .T.C
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
          cursteps = cursteps + 1
        end do
        curministeps = curministeps + steps
        if (present(upgraded)) then
          steps = 0
          cursteps = cursteps+1
          do cr = 1, 2
            dops(cr) = floor(min(sqrt(CC(cr)%n-k*1.0d0), CC(cr)%n/(k*1.0d0)))
            call p(cr)%pert(CC(cr)%n)
            call CC0(cr)%sort(p(cr))
            p(cr) = p(cr)%reverse()
          end do
          i20 = 0
          j20 = 0
          i21 = 0
          j21 = 0
          CIJ = 1.0d0
          do cr = 1, dops(1)
            j2 = floor(p(1)%d(cr) + 0.5)
            CJ = CC(1)%submatrix(j2,k,j2,1)
            call CJ%maxelement(i1, i2)
            call this%swap(1, i2, j2)
            call CJ%deinit()
            call CJ%init(k, dops(2) + k)
            do i1 = 1, k
              do j1 = 1, k
                CJ%d(i1,j1) = A%d(i1,j1)
              end do
            end do
            do j1 = k+1, k + dops(2)
              j2 = floor(p(2)%d(j1-k)+0.5)
              do i1 = 1, k
                CJ%d(i1, j1) = this%d(i1, j2)
              end do
            end do
            CJ = A .Id. CJ
            CJ = .T.CJ
            call CJ%maxelement(i1,j1)
            if (abs(CJ%d(i1,j1))*CC0(1)%d(cr) > CIJ) then
              i20 = i2
              i21 = floor(p(1)%d(cr) + 0.5)
              j20 = j1
              j21 = floor(p(2)%d(i1 - k) + 0.5)
              !Почему остаются 1 при сортировке?
              CIJ = abs(CJ%d(i1,j1))*CC0(1)%d(cr)
            end if
            j2 = floor(p(1)%d(cr) + 0.5)
            call this%swap(1, i2, j2)
            call CJ%deinit()
          end do
          if (CIJ > 1.0d0) then
            print *, CIJ, i20, i21, j20, j21
            steps = 1
            call this%swap(1, i20, i21)
            call this%swap(2, j20, j21)
            if ((present(pert1)) .and. (present(pert2))) then
              tmp = pert1%d(i20)
              pert1%d(i20) = pert1%d(i21)
              pert1%d(i21) = tmp
              tmp = pert2%d(j20)
              pert2%d(j20) = pert2%d(j21)
              pert2%d(j21) = tmp
            end if
          end if
          do cr = 1, 2
            call p(cr)%deinit()
          end do
        end if
        if (present(maxsteps)) then
          if (cursteps >= maxsteps) then
            exit
          end if
        end if
      end do
      if (present(bigsteps)) then
        bigsteps = cursteps
      end if
      if (present(totsteps)) then
        totsteps = curministeps
      end if
    end
    
    subroutine mtrx_maxvol25(this, k, l, pert1, pert2, bigsteps)
      Class(Mtrx) :: this
      Type(Vector), optional :: pert1, pert2
      Integer(4), intent(in) :: k, l
      Integer(4), intent(out), optional :: bigsteps
      Integer(4) :: cursteps, curministeps, steps, prevsteps
      Integer(4) n, m, cr, maxst, maxpert, curpert
      if (k .eq. l) then
        if (present(bigsteps)) then
          call this%maxvol(k, pert1, pert2, bigsteps)
        else
          call this%maxvol(k, pert1, pert2)
        end if
        return
      end if
      n = this%n
      m = this%m
      cursteps = 1
      maxst = 20
      maxpert = 2*max(k,l)
      if ((k > n) .or. (l > m)) then
        print *, "error in maxvol25"
      end if
      !Почему-то хуже
      !if (k > l) then
      !  call this%premaxvol(2,k,l,pert2)
      !else
      !  call this%premaxvol(1,l,k,pert1)
      !end if
      steps = 1
      prevsteps = 1
      do while (steps > 0)
        steps = 0
        do cr = 1, 2
          curpert = 0
          if (cr .eq. 1) then
            if (k > l) then
              call this%maxvol251(cr,l,k,0,pert1,pert2,curpert)
            else
              call this%dominantr(cr,k,l,pert1,curpert)
            end if
          else
            if (k > l) then
              call this%dominantr(cr,l,k,pert2,curpert)
            else
              call this%maxvol251(cr,k,l,0,pert1,pert2,curpert)
            end if
          end if
          
          steps = steps + curpert
          prevsteps = prevsteps + curpert
          if ((steps == 0) .and. (cr == 2)) then
            exit
          end if
          if ((prevsteps == 0) .and. (cr == 1)) then
            steps = 0
            exit
          end if
          if (cursteps >= maxst) then
            steps = 0
            exit
          end if
          if (cr == 1) then
            prevsteps = 0
          end if
          cursteps = cursteps + 1
        end do
        curministeps = curministeps + steps
      end do
      if (present(bigsteps)) then
        bigsteps = cursteps
      end if
    end
    
    subroutine mtrx_hmaxvol2(this, tin, k, l, pert, CIN, addmask)
      Class(Mtrx) :: this
      Type(Mtrx) :: C
      Type(Vector), optional :: pert
      Integer(4), optional :: CIN
      Type(Mtrx), intent(in), optional :: addmask
      Type(Vector) :: LB, CI, CJ
      Integer(4), intent(in) :: k, l
      Integer(4), intent(in) :: tin
      Integer(4) t
      Integer(4) n, i1, cr
      Integer(4) ij(1)
      DOUBLE PRECISION tau, ls, alpha
      DOUBLE PRECISION, allocatable :: work(:)
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
      Allocate(work(n-k))
      LB = (C .dot. C)*evec(k)
      do cr = 1, l-k
        ij = maxloc(LB%d)
        i1 = ij(1)
        call this%swap(t, i1+k, cr+k)
        if (present(pert)) then
          call pert%swap(i1+k, cr+k)
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
    
    subroutine mtrx_maxvol2(this, t, k, l, pert, addmask)
      Class(Mtrx) :: this
      Type(Mtrx) :: A, C, CI, CB, CBI, CS
      Type(Vector), optional :: pert
      Type(Mtrx), intent(in), optional :: addmask
      Type(Vector) :: LB
      Integer(4), intent(in) :: k, l
      Integer(4), intent(in) :: t
      Integer(4) n, m, i, j, i1, j1, cr, tm
      Integer(4) ij(1)
      DOUBLE PRECISION tmp, ls
      n = this%n
      m = this%m
      if (((l > n) .and. (t == 1)) .or. ((l > m) .and. (t .ne. 1))) then
        print *, "error in maxvol2"
      end if
      if (l < k) then
        print *, "error in maxvol2"
      end if
      call LB%init(max(m,n))
      A = this%submatrix(k, k)
      if (t .eq. 1) then
      !Добавляем строки
        tm = 1
        C = this%submatrix(n, k)
        C = C .dI. A
      else
      !Добавляем столбцы
        n = m
        tm = 2
        C = this%submatrix(k, n)
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
        if (present(pert)) then
          if (tm == 1) then
            tmp = pert%d(i)
            pert%d(i) = pert%d(cr)
            pert%d(cr) = tmp
          else
            tmp = pert%d(i)
            pert%d(i) = pert%d(cr)
            pert%d(cr) = tmp
          end if
        end if
        if (present(addmask)) then
          call addmask%swap(tm, i, cr)
        end if
        ls = LB%d(i)
        CI = CB%submatrix(i, l, i, 1)
        CBI = (.T.CI)/ls
        CS = CB * CBI
        !CB = CB - (CS * CI)
        call CB%update1(-1.0d0, CS, CI)
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
    
!     subroutine mtrx_dominantc(this, t, k, l, nai, pert1, pert2, steps, maxstepsin)
!       Class(Mtrx) :: this
!       Type(Mtrx) :: A, C, CI, CJ, CJ1, CJ2, CI1, CI2, CB, CBI, CBJ, CBJ2
!       Type(Vector), optional :: pert1, pert2
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
!       !Изменяем строки
!         tm = 1
!         A = this%submatrix(l, k)
!         C = this%submatrix(n, k)
!         C = C .dI. A
!       else
!       !Изменяем столбцы
!         n = m
!         tm = 2
!         A = this%submatrix(k, l)
!         C = this%submatrix(k, n)
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
!         if ((present(pert1)) .and. (present(pert2))) then
!           if (tm == 1) then
!             tmp = pert1%d(i)
!             pert1%d(i) = pert1%d(j)
!             pert1%d(j) = tmp
!           else
!             tmp = pert2%d(i)
!             pert2%d(i) = pert2%d(j)
!             pert2%d(j) = tmp
!           end if
!         end if
!         k11 = CB%d(i,j)/ro
!         k12 = (1.0d0 - CB%d(j,j))/ro
!         k21 = LB%d(i)/ro
!         k22 = k11
!         CJ = CB%submatrix(j, l, j, 1)
!         CJ%d(1, j) = CJ%d(1,j)-1.0d0
!         CJ1 = CJ*k21
!         CJ2 = CJ*k22
!         CI = CB%submatrix(i, l, i, 1)
!         CBI = CB*(.T.CI)
!         CI%d(1, j) = CI%d(1,j)-1.0d0
!         CI1 = CI*k11
!         CI2 = CI*k12
!         CBJ = CB%submatrix(n, j, 1, j)
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

subroutine mtrx_dominantc(this, t, k, l, nai, pert1, pert2, steps, maxstepsin)
      Class(Mtrx) :: this
      Type(Mtrx) :: A, CI, CI1, CB, CBI, B
      Type(Vector), optional :: pert1, pert2
      Type(Vector) :: LB, LC
      Integer(4), intent(in), optional :: nai
      Integer(4), intent(in) :: k, l
      Integer(4), intent(in) :: t
      Integer(4), intent(out), optional :: steps
      Integer(4), intent(in), optional :: maxstepsin
      Integer(4) n, m, i, j, i1, j1, tm, maxsteps, cursteps, na
      Integer(4) ij(2)
      DOUBLE PRECISION tmp, k11, k12
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
      !Изменяем строки
        tm = 1
        A = this%submatrix(l, k)
        CB = this%submatrix(n, k)
        CB = CB .dI. A
      else
      !Изменяем столбцы
        n = m
        tm = 2
        A = this%submatrix(k, l)
        CB = this%submatrix(k, n)
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
        if ((present(pert1)) .and. (present(pert2))) then
          if (tm == 1) then
            tmp = pert1%d(i)
            pert1%d(i) = pert1%d(j)
            pert1%d(j) = tmp
          else
            tmp = pert2%d(i)
            pert2%d(i) = pert2%d(j)
            pert2%d(j) = tmp
          end if
        end if
        
        k11 = LB%d(i)
        CBI = CB%submatrix(i, l, i, 1)/k11
        CI = CB*(.T.CBI)
        do i1 = 1, n
          LB%d(i1) = LB%d(i1) - k11*CI%d(i1,1)**2
        end do
        !CB = CB - CI*CB%submatrix(i, l, i, 1)
        call CB%update1(-1.0d0,CI,CB%submatrix(i, l, i, 1))
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
        !CB = CB + k12*CI*CB%submatrix(i, l, i, 1)
        call CB%update1(k12,CI,CB%submatrix(i, l, i, 1))
        !CB = CB + k12*CI*CB%submatrix(j, l, j, 1)
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
    
    subroutine mtrx_maxvol251(this, t, k, l, nai, pert1, pert2, steps, maxstepsin)
      Class(Mtrx) :: this
      Type(Mtrx) :: A, CB, B
      Type(Vector) :: CI, CJ, CJ1, CJ2, CI1, CI2, CBI, CBJ
      Type(Vector), optional :: pert1, pert2
      Type(Vector) :: LB, LC
      Integer(4), intent(in), optional :: nai
      Integer(4), intent(in) :: k, l
      Integer(4), intent(in) :: t
      Integer(4), intent(out), optional :: steps
      Integer(4), intent(in), optional :: maxstepsin
      Integer(4) n, m, i, j, i1, j1, tm, maxsteps, cursteps, na
      Integer(4) ij(2)
      DOUBLE PRECISION tmp, k11, k12, k21, k22, li
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
      !Изменяем строки
        tm = 1
        A = this%submatrix(l, k)
        CB = this%submatrix(n, k)
        CB = CB .dI. A
      else
      !Изменяем столбцы
        n = m
        tm = 2
        A = this%submatrix(k, l)
        CB = this%submatrix(k, n)
        CB = A .Id. CB
        CB = .T.CB
      end if
      LB = (CB .dot. CB) * evec(l)
      do i1 = l+1, n
        LB%d(i1) = LB%d(i1) + 1.0d0
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
      maxsteps = 2*l*2
      if (present(maxstepsin)) then
        maxsteps = maxstepsin
      end if
      cursteps = 0
      do while (B%d(i,j) > 1.0d0)
        if (cursteps >= maxsteps) exit
        call this%swap(tm, i, j)
        if ((present(pert1)) .and. (present(pert2))) then
          if (tm == 1) then
            tmp = pert1%d(i)
            pert1%d(i) = pert1%d(j)
            pert1%d(j) = tmp
          else
            tmp = pert2%d(i)
            pert2%d(i) = pert2%d(j)
            pert2%d(j) = tmp
          end if
        end if
        k11 = CB%d(i,j)/B%d(i,j)
        k12 = (1.0d0 - CB%d(j,j))/B%d(i,j)
        k21 = LB%d(i)/B%d(i,j)
        k22 = k11
        CJ = tovec(CB%submatrix(j, l, j, 1))
        CJ%d(j) = CJ%d(j)-1.0d0
        CJ1 = CJ*k21
        CJ2 = CJ*k22
        CI = tovec(CB%submatrix(i, l, i, 1))
        CBI = CB*CI
        CI%d(j) = CI%d(j)-1.0d0
        CI1 = CI*k11
        CI2 = CI*k12
        CBJ = tovec(CB%submatrix(n, j, 1, j))
        li = LB%d(i)
        call CB%update1v(1.0d0,CBJ,CJ1-CI1)
        call CB%update1v(-1.0d0,CBI,CJ2+CI2)
        call CB%swap(1, i, j)
        call CBJ%update1(-k11/k21, CBI)
        call LB%update1(-1.0d0/li, CBI .dot. CBI)
        call LB%update1(k21, CBJ .dot. CBJ)
        call LB%swap(i,j)
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
    
    subroutine mtrx_premaxvol(this, t, l, maxr, pert)
      Class(Mtrx) :: this
      Integer(4), intent(in) :: t, l
      Type(Vector), intent(inout), optional :: pert
      Integer(4), intent(in), optional :: maxr
      Type(Vector) c, w, piv
      Type(Mtrx) mat, mat2, da, AB, q, A, dapr2, bc
      Integer(4) j, i, rank
      Integer(4) n, m, maxrank
      Integer(4) ij(1)
      DOUBLE PRECISION tmp, aa
      
      if (((l > this%m) .and. (t == 1)) .or. ((l > this%n) .and. (t .ne. 1))) then
        print *, "error in premaxvol"
      end if
      if (l < maxr) then
        print *, "error in premaxvol"
      end if
      
      if (t .eq. 1) then
        da = .T.(this%submatrix(this%n, l))
      else
        da = this%submatrix(l, this%m)
      end if
      
      n = da%n
      m = da%m
      
      call c%init(m)
      call w%init(maxr)
      call AB%init(maxr,m)
      call A%init(maxr,maxr)
      call piv%init(m)
      do j = 1, m
        mat = da%submatrix(n, j, 1, j)
        mat = (.T.mat) * mat
        c%d(j) = mat%d(1,1)
        piv%d(j) = j
      end do
      if (present(pert)) then
        piv = pert
      end if
      rank = 0
      if (present(maxr)) then
        maxrank = maxr
      else
        maxrank = l
      end if
      ij = maxloc(c%d)
      j = ij(1)
      
      q = eye(n)
      
      do while (rank < maxrank)
      
        rank = rank + 1
        
        !Перемещаем j на место rank
        tmp = piv%d(j)
        piv%d(j) = piv%d(rank)
        piv%d(rank) = tmp
        tmp = c%d(j)
        c%d(j) = c%d(rank)
        c%d(rank) = tmp
        call AB%swap(2, j, rank)
        call this%swap(t, j, rank)
        call da%swap(2, j, rank)
        aa = sqrt(c%d(rank))
        if (t .eq. 1) then
          bc = .T.(this%submatrix(rank, n, rank, 1))
        else
          bc = this%submatrix(n, rank, 1, rank)
        end if
        if (rank > 1) then
          mat = A%submatrix(rank-1,rank-1)*AB%submatrix(rank-1,rank,1,rank)
        end if
        do i = 1, rank-1
          A%d(i, rank-1) = mat%d(i,1)
        end do
        A%d(rank,rank) = aa
        if (rank > 1) then
          mat2 = (bc - q%submatrix(n,rank-1)*mat)/aa
        else
          mat2 = bc/aa
        end if
        do j = 1, n
          q%d(j,rank) = mat2%d(j,1)
        end do
        if (t .eq. 1) then
          mat = .T.(this%submatrix(m, n, rank, 1))
        else
          mat = this%submatrix(n, m, 1, rank)
        end if
        if (rank > 1) then
          mat2 = q%submatrix(n,rank-1)
          mat2 = ((.T.bc) - ((.T.bc)*mat2)*(.T.mat2))/aa
        else
          mat2 = (.T.bc)/aa
        end if
        dapr2 = mat2 * mat
        
        !Пересчет c
        c%d(rank) = 0
        do i = rank + 1, m
          c%d(i) = c%d(i) - dapr2%d(1,i-rank+1)**2
        end do
        
        !Пересчет w
        w%d(rank) = 1.0d0/aa**2
        do i = 1, rank - 1
          w%d(i) = w%d(i) + AB%d(i, rank)**2*w%d(rank)
        end do
        
        !Пересчет AB
        do i = rank+1, m
          AB%d(rank, i) = dapr2%d(1,i-rank)/aa
        end do
        if (rank > 1) then
          mat = AB%submatrix(rank-1,rank,1,rank)*(dapr2/aa)
        end if
        do i = rank+1, m
          do j = 1, rank-1
            AB%d(j, i) = AB%d(j,i) - mat%d(j, i-rank)
          end do
        end do
        if (rank > 1) then
          call mat%deinit()
        end if

        ij = maxloc(c%d)
        j = ij(1)
      end do
      if (present(pert)) then
        pert = piv
      end if
    end
    
    subroutine mtrx_maxvol252(this, t, rank, l, pert, steps, maxministepsin)
      Class(Mtrx) :: this
      Type(Vector), intent(inout), optional :: pert
      Integer(4), intent(in) :: t, rank, l
      Integer(4), intent(out), optional :: steps
      Integer(4), intent(in), optional :: maxministepsin
      Type(Vector) c, w, piv
      Type(Mtrx) mat, mat2, da, dapr1, dapr2, dapr3, AB, u1, u2, u, A, bnew, cnew, q, bc
      Integer(4) j, i, i1, j1
      Integer(4) n, m, maxministeps, curministeps
      Integer(4) ij(2)
      DOUBLE PRECISION tmp, beta, ro, alpha, mu, f, dapr0, aa, bnew0, cnew0
      
      if (rank .eq. l) then
        !Пока такого нет:
        !call this%maxvol(t, rank)
        return
      end if
      
      if (t .eq. 1) then
        da = .T.(this%submatrix(this%n, l))
      else
        da = this%submatrix(l, this%m)
      end if
      
      n = da%n
      m = da%m
      
      f = 1.0d0
      call c%init(m)
      call w%init(min(m,n))
      call AB%init(min(m,n),m)
      call piv%pert(m)
      if (present(pert)) then
        piv = pert
      end if
      maxministeps = floor(log(1.0d0*m)/log(2.0d0))+2*l
      if (present(maxministepsin)) then
        maxministeps = maxministepsin
      end if
      
      mat = da%submatrix(n, rank)
      call mat%qr(u, A, 1)
      q = u%submatrix(n,rank+1)
      do j = 1, m-rank
        mat = da%submatrix(n, j+rank, 1, j+rank)
        mat = (.T.mat) * mat
        c%d(j+rank) = mat%d(1,1)
      end do
      da = (.T.u)*da
      do j = 1, m-rank
        mat = da%submatrix(rank, j+rank, 1, j+rank)
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
      AB = mat2*da%submatrix(rank,m)
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
          
          !Меняем местами i1 и rank
          do i = i1+1, rank
            tmp = w%d(i-1)
            w%d(i-1) = w%d(i)
            w%d(i) = tmp
            tmp = piv%d(i-1)
            piv%d(i-1) = piv%d(i)
            piv%d(i) = tmp
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
            mat2 = mat * A%submatrix(i+1,rank,i,i)
            do j = i, rank
              A%d(i,j) = mat2%d(1,j-i+1)
              A%d(i+1,j) = mat2%d(2,j-i+1)
            end do
            mat2 = q%submatrix(n,i+1,1,i) * (.T.mat)
            do j = 1, n
              q%d(j,i) = mat2%d(j,1)
              q%d(j,i+1) = mat2%d(j,2)
            end do
            call mat2%deinit()
            call mat%deinit()
          end do
          
          !Меняем местами: rank, rank+1, j1 => j1, rank, rank+1
          tmp = piv%d(j1)
          piv%d(j1) = piv%d(rank+1)
          piv%d(rank+1) = tmp
          call this%swap(t, j1, rank+1)
          tmp = c%d(j1)
          c%d(j1) = c%d(rank+1)
          c%d(rank+1) = tmp
          aa = sqrt(c%d(rank+1))
          if (t .eq. 1) then
            bc = .T.(this%submatrix(rank+1, n, rank+1, 1))
          else
            bc = this%submatrix(n, rank+1, 1, rank+1)
          end if
          call AB%swap(2, j1, rank+1)
          dapr1 = A%d(rank,rank)*AB%submatrix(rank,m,rank,rank+2)
          dapr0 = A%d(rank,rank)
          if (rank > 1) then
            dapr3 = A%submatrix(rank-1,rank,1,rank)
          end if
          mat = A*AB%submatrix(rank,rank+1,1,rank+1)
          do i = 1, rank
            A%d(i, rank) = mat%d(i,1)
          end do
          mat2 = (bc - q%submatrix(n,rank)*mat)/aa
          do j = 1, n
            q%d(j,rank+1) = mat2%d(j,1)
          end do
          if (t .eq. 1) then
            mat = .T.(this%submatrix(m, n, rank+1, 1))
          else
            mat = this%submatrix(n, m, 1, rank+1)
          end if
          mat2 = q%submatrix(n,rank)
          mat2 = ((.T.bc) - ((.T.bc)*mat2)*(.T.mat2))/aa
          dapr2 = mat2 * mat
          call mat%deinit()
          call mat2%deinit()
          tmp = piv%d(rank)
          piv%d(rank) = piv%d(rank + 1)
          piv%d(rank + 1) = tmp
          call this%swap(t, rank, rank + 1)
          
          !Обнуляем столбец rank под диагональю (одно лишнее число)
          call givens(alpha, beta, A%d(rank,rank), aa)
          call mat%init(2,2)
          mat%d(1,1) = alpha
          mat%d(1,2) = -beta
          mat%d(2,1) = beta
          mat%d(2,2) = alpha
          mat2 = q%submatrix(n,rank+1,1,rank) * (.T.mat)
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
          mat2 = mat%submatrix(2,1) * u + mat%submatrix(2,2,1,2) * u2
          cnew0 = mat2%d(2,1)
          cnew = mat2%submatrix(2,m-rank,2,2)
          bnew0 = mat2%d(1,1)
          bnew = mat2%submatrix(1,m-rank,1,2)
          A%d(rank,rank) = sqrt(aa**2 + A%d(rank,rank)**2)
          call mat2%deinit()
          call u%deinit()
          call u2%deinit()
          call mat%deinit()
          
          !Пересчет c
          c%d(rank+1) = cnew0**2
          do i = rank + 2, m
            c%d(i) = c%d(i) + cnew%d(1,i-rank-1)**2 - dapr2%d(1,i-rank)**2
          end do
          
          !Пересчет w
          if (rank > 1) then
            mat = A%submatrix(rank-1,rank-1)
            u = mat%rtsolve(dapr3)
          end if
          w%d(rank) = 1.0d0/A%d(rank,rank)**2
          if (rank > 1) then
            u2 = AB%submatrix(rank-1,rank+1,1,rank+1) + AB%d(rank,rank+1)*u
          end if
          do i = 1, rank - 1
            w%d(i) = w%d(i) + w%d(rank)*u2%d(i,1)**2 - u%d(i,1)**2/dapr0**2
          end do
          
          !Пересчет AB
          mu = AB%d(rank,rank+1)
          if (rank > 1) then
            u1 = AB%submatrix(rank-1,rank+1,1,rank+1)
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
      if (present(pert)) then
        pert = piv
      end if
      if (present(steps)) then
        steps = curministeps
      end if
    end
    
    subroutine mtrx_dominantr(this, t, r, l, pert, steps, maxstepsin, AIout, ABout)
      Class(Mtrx) :: this
      Type(Vector), intent(inout), optional :: pert
      Integer(4), intent(in) :: t, r, l
      Integer(4), intent(out), optional :: steps
      Integer(4), intent(in), optional :: maxstepsin
      Type(Mtrx), intent(out), optional :: AIout, ABout
      Type(Mtrx) da, db
      Type(Vector) c, w
      Type(Mtrx) AB, AI, B
      Type(Mtrx) mat
      Type(Mtrx) q
      Integer(4) info
      Integer(4) ij(2), i1, j1, i
      DOUBLE PRECISION ro, l1, abij
      Integer(4) cursteps, maxsteps
      Type(Vector) c1, c2, abj, aig, aii, abi
      Type(Vector) v
      
      if (t .eq. 1) then
        da = .T.(this%subarray(r, l))
        db = .T.(this%subarray(this%n, l, r+1, 1))
      else
        da = this%subarray(l, r)
        db = this%subarray(l, this%m, 1, r+1)
      end if
      call da%qr(q, AI)
      B = (.T.q) * db
      AB = AI%rtsolve(B)
      call dtrtri('U', 'N', r, AI%d, r, info)
      w = (AI .dot. AI) * evec(r)
      c = (evec(l) * (db .dot. db)) - (evec(r) * (B .dot. B))

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
      do while ((abs(ro) > 1.0d0) .and. (cursteps < maxsteps))
        !Введение переменных
        abij = AB%d(i1,j1)
        l1 = sqrt(c%d(j1))
        aii = tovec(AI%subarray(i1, r, i1, 1))
        aig = AI*aii
        abi = tovec(AB%subarray(i1, AB%m, i1, 1))
        abj = tovec(AB%subarray(r, j1, 1, j1))
        c1 = ((tovec(db%subarray(l, j1, 1, j1)) - da*abj) / l1) * db
        !swap columns
        call c2%copy(c1)
        c2%d(j1) = 0.0d0
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
        if (present(pert)) then
          call pert%swap(i1,j1+r)
        end if
        !Пересчет c
        c2 = c2 * (-abij/ro)
        call c2%update1(l1/ro, abi)
        c = c - (c1 .dot. c1) + (c2 .dot. c2)
        c1%d(j1) = 0.0d0
        !Перечет AI
        call v%copy(abj)
        call v%update1((ro - abij) / w%d(i1), aig)
        v%d(i1) = v%d(i1) - 1.0d0
        call AI%update1v(-1.0d0/(v%d(i1) + 1.0d0), v, aii)
        !Пересчет w
        w = w + (v .dot. (v * (w%d(i1) / (v%d(i1) + 1.0d0)**2) - aig * (2.0d0 / (v%d(i1) + 1.0d0))))
        !Пересчет AB
        call AB%update1v(-1.0d0/(v%d(i1) + 1.0d0), v, abi)
        aig = aig * (l1/(ro + abij))
        call aig%update1(-aig%d(i1)/(v%d(i1) + 1.0d0), v)
        call AB%update1v(1.0d0, aig, c1-c2)

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
        !Добавляем строки
          A = this%submatrix(cr-1, l)
          C = this%submatrix(n, l)
          C = C .dI. A
        else
        !Добавляем столбцы
          A = this%submatrix(l, cr-1)
          C = this%submatrix(l, n)
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
    
    !Решить, оставить это или maxvol25
    !Пока pert не сохраняются
    subroutine mtrx_maxvolrect(this, k, l, pert1in, pert2in, bigsteps, totsteps)
      Class(Mtrx) :: this
      Integer(4), intent(in) :: k, l
      Type(Vector) :: pert1in, pert2in
      Integer(4), optional :: bigsteps, totsteps
      Type(Vector) pert1, pert2
      Integer(4) cursteps, curministeps, steps, prevsteps
      Integer(4) cr, maxst, maxpert, curpert
      Type(Mtrx) cols, rows
      Integer(4) m, n
      n = this%n
      m = this%m
      if (k .eq. l) then
        cursteps = bigsteps
        call this%maxvol(k, pert1, pert2, bigsteps, totsteps, cursteps)
        return
      end if
      call cols%init(this%n, l)
      call rows%init(k, this%m)
      cursteps = 1
      maxst = 20
      if ((present(bigsteps)) .and. (bigsteps > 0)) then
        maxst = bigsteps
      end if
      maxpert = 2*max(k,l)
      if ((present(totsteps)) .and. (totsteps > 0)) then
        maxpert = totsteps
      end if
      if ((k > n) .or. (l > m)) then
        print *, "error in maxvol-rect"
      end if
      steps = 1
      prevsteps = 1
      do while (steps > 0)
        steps = 0
        do cr = 1, 2
          curpert = 0
          if (cr .eq. 1) then
            call pert1%pert(this%n)
            call pert2%pert(this%m)
            cols = this%submatrix(this%n, l)
            if (k > l) then
              call cols%maxvol251(1,l,k,0,pert1,pert2,curpert,k)
            else
              call cols%dominantr(1,k,l,pert1,curpert)
            end if
            call this%pertrows(pert1, 1, k)
            call this%pertcols(pert2, 1, l)
            call pert1in%pertapp(pert1, k)
            call pert2in%pertapp(pert2, l)
            call pert1%deinit()
            call pert2%deinit()
          else
            call pert1%pert(this%n)
            call pert2%pert(this%m)
            rows = this%submatrix(k, this%m)
            if (k > l) then
              call rows%dominantr(2,l,k,pert2,curpert)
            else
              call rows%maxvol251(2,k,l,0,pert1,pert2,curpert,l)
            end if
            call this%pertrows(pert1, 1, k)
            call this%pertcols(pert2, 1, l)
            call pert1in%pertapp(pert1, k)
            call pert2in%pertapp(pert2, l)
            call pert1%deinit()
            call pert2%deinit()
          end if
          steps = steps + curpert
          prevsteps = prevsteps + curpert
          if ((steps == 0) .and. (cr == 2)) then
            exit
          end if
          if ((prevsteps == 0) .and. (cr == 1)) then
            steps = 0
            exit
          end if
          if (cursteps >= maxst) then
            steps = 0
            exit
          end if
          if (cr == 1) then
            prevsteps = 0
          end if
          cursteps = cursteps + 1
        end do
        curministeps = curministeps + steps
      end do
      if (present(bigsteps)) then
        bigsteps = cursteps
      end if
      if (present(totsteps)) then
        totsteps = curministeps
      end if
    end
    
    !Пока не сохраняет финальные pert
    subroutine mtrx_maxvolproj(this, r, k, l, pert1, pert2, bigsteps, totsteps)
      Class(Mtrx) :: this
      Integer(4), intent(in) :: r, k, l
      Type(Vector), intent(out) :: pert1, pert2
      Integer(4), optional :: bigsteps, totsteps
      Integer(4) :: cursteps, curministeps, steps, prevsteps
      Integer(4) cr, maxst, maxpert, curpert
      Integer(4) m, n
      Type(Mtrx) cols, rows
      Type(Mtrx) ai, u, s, v
      n = this%n
      m = this%m
      if ((k .eq. r) .or. (l .eq. r)) then
        call this%maxvolrect(k, l, pert1, pert2, bigsteps, totsteps)
        return
      end if
      call pert1%pert(this%n)
      call pert2%pert(this%m)
      call this%maxvol(r, pert1, pert2, cursteps, curministeps, 4)
      call pert1%deinit()
      call pert2%deinit()
      call pert1%pert(this%n)
      call pert2%pert(this%m)
      rows = this%submatrix(r, this%m)
      call rows%maxvol251(2, r, l, 0, pert1, pert2, curpert, l)
      call this%pertcols(pert2, 1, l)
      call pert1%deinit()
      call pert2%deinit()
      call pert1%pert(this%n)
      call pert2%pert(this%m)
      cols = this%submatrix(this%n, l) * (.I.this%submatrix(r, l))
      call cols%maxvol251(1, r, k, 0, pert1, pert2, curpert, k)
      call this%pertrows(pert1, 1, k)
      call pert1%deinit()
      call pert2%deinit()
      call cols%deinit()
      call rows%deinit()
      call cols%init(this%n, l)
      call rows%init(k, this%m)
      cursteps = 1
      maxst = 4
      if ((present(bigsteps)) .and. (bigsteps > 0)) then
        maxst = bigsteps
      end if
      maxpert = max(k,l)
      if ((present(totsteps)) .and. (totsteps > 0)) then
        maxpert = totsteps
      end if
      if ((k > n) .or. (l > m)) then
        print *, "error in maxvolproj"
      end if
      steps = 1
      prevsteps = 1
      do while (steps > 0)
        steps = 0
        do cr = 1, 2
          curpert = 0
          if (cr .eq. 2) then
            cols = this%submatrix(this%n, l)
            ai = this%submatrix(k,l)
            call ai%svd(u,s,v)
            v = v%submatrix(r, l)
            cols = cols*(.T.v)
            call pert1%pert(this%n)
            call pert2%pert(this%m)
            call cols%maxvol251(1,r,k,0,pert1,pert2,curpert,k)
            call this%pertrows(pert1, 1, k)
            call pert1%deinit()
            call pert2%deinit()
          else
            rows = this%submatrix(k, this%m)
            ai = this%submatrix(k,l)
            call ai%svd(u,s,v)
            u = u%submatrix(k, r)
            rows = (.T.u)*rows
            call pert1%pert(this%n)
            call pert2%pert(this%m)
            call rows%maxvol251(2,r,l,0,pert1,pert2,curpert,l)
            call this%pertcols(pert1, 1, l)
            call pert1%deinit()
            call pert2%deinit()
          end if
          steps = steps + curpert
          prevsteps = prevsteps + curpert
          if ((steps == 0) .and. (cr == 2)) then
            exit
          end if
          if ((prevsteps == 0) .and. (cr == 1)) then
            steps = 0
            exit
          end if
          if (cursteps >= maxst) then
            steps = 0
            exit
          end if
          if (cr == 1) then
            prevsteps = 0
          end if
          cursteps = cursteps + 1
        end do
        curministeps = curministeps + steps
      end do
      if (present(bigsteps)) then
        bigsteps = cursteps
      end if
      if (present(totsteps)) then
        totsteps = curministeps
      end if
    end
  
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
    
    subroutine mtrx_pertrows(this, pert, t, kin)
      Class(Mtrx) :: this
      Type(Mtrx) mat
      Type(Vector), intent(in) :: pert
      Integer(4) n, m, i, j, tmp
      Integer(4), intent(in) :: t
      Integer(4), intent(in), optional :: kin
      Integer(4) k
      Real(8) tmpr
      n = this%n
      m = this%m
      k = n
      if (present(kin)) then
        k = kin
      end if
      if (pert%n .ne. n) then
        print *, "error in pert_rows_mtrx", pert%n, n
        return
      end if
      call mat%init(n,m)
      do i = 1, k
        tmp = floor(pert%d(i) + 0.5)
        if (t == 1) then
          do j = 1, m
            tmpr = mat%d(i,j)
            mat%d(i,j) = this%d(tmp,j)
            if (tmp > k) then
              this%d(tmp,j) = tmpr
            end if
          end do
        else
          do j = 1, m
            mat%d(tmp,j) = this%d(i,j)
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
    
    subroutine mtrx_pertcols(this, pert, t, kin)
      Class(Mtrx) :: this
      Type(Mtrx) mat
      Type(Vector), intent(in) :: pert
      Integer(4) n, m, i, j, tmp, k
      Integer(4), intent(in) :: t
      Integer(4), intent(in), optional :: kin
      Real(8) tmpr
      n = this%n
      m = this%m
      k = m
      if (present(kin)) then
        k = kin
      end if
      if (pert%n .ne. m) then
        print *, "error in pert_cols_mtrx", pert%n, m
        return
      end if
      call mat%init(n,m)
      do i = 1, k
        tmp = floor(pert%d(i) + 0.5)
        if (t == 1) then
          do j = 1, n
            tmpr = mat%d(j,i)
            mat%d(j,i) = this%d(j,tmp)
            if (tmp > k) then
              this%d(j, tmp) = tmpr
            end if
          end do
        else
          do j = 1, n
            mat%d(j,tmp) = this%d(j,i)
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
    
    function mtrx_vol(this, rank) Result(res)
      Class(Mtrx), intent(in) :: this
      Integer(4), intent(in), optional :: rank
      Type(Mtrx) :: u, s, v
      DOUBLE PRECISION res
      Integer(4) n, i
      if (.not. present(rank)) then
        if (this%n > this%m) then
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
      res%d = this%d(n1:n,m1:m)
    end
  
    function mtrx_subm(this, n, m, k, l) Result(res)
    !k - начальная строка; l - начальный столбец. n и m - конечные
      Class(Mtrx), intent(in) :: this
      Integer(4), intent(in) :: n, m
      Integer(4), intent(in), optional :: k, l
      Integer(4) n1, m1
      Type(Mtrx) :: res
      Integer(4) i, j
      if ((.not. present(k)) .or. (.not. present(l))) then
        n1 = 1
        m1 = 1
      else
        n1 = k
        m1 = l
      end if
      if ((n1 > n) .or. (m1 > m)) then
        print *, "error submatrix_mtrx too big", n, m, n1, m1
        !call backtrace()
        stop
      end if
      call res%init(n - n1 + 1, m - m1 + 1)
      if ((n <= this%n) .and. (m <= this%m)) then
        do j = m1, m
          do i = n1, n
            res%d(i - n1 + 1, j - m1 + 1) = this%d(i, j)
          end do
        end do
      else
        print *, "error submatrix_mtrx", n, m, this%n, this%m
        !call backtrace()
        stop
      end if
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
      DOUBLE PRECISION, Allocatable :: ds(:), du(:,:), dvt(:,:), da(:,:)
      DOUBLE PRECISION, Allocatable :: work(:)
      Integer(4) Lwork, info
      Integer(4) n, m, i, j
      n = this%n
      m = this%m
      Lwork = 5 * (m + n)
      Allocate(work(Lwork))
      Allocate(da(n, m))
      Allocate(ds(n))
      Allocate(du(n, n))
      Allocate(dvt(m, m))
      do j = 1, m
        do i = 1, n
          da(i, j) = this%d(i, j)
        end do
      end do
      call dgesvd('A', 'A', n, m, da, n, ds, du, n, dvt, m, work, Lwork, info)
      Deallocate(work)
      Deallocate(da)
      call u%init(n, n)
      call s%init(n, m)
      call vt%init(m, m)
      call u%set(du)
      call vt%set(dvt)
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
!       Type(Mtrx) mat
!       Integer(4) n, m
!       n = this%n
!       m = this%m
!       if ((n .ne. b%n) .or. (b%m .ne. 1) .or. (n .ne. m)) then
!         print *, 'error rtsolve', n, m, b%n, b%m
!       end if
!       call res%init(n,1)
!       res%d(n,1) = b%d(n,1)/this%d(n,n)
!       do i = n-1, 1, -1
!         mat = this%submatrix(i,n,i,i+1)*res%submatrix(n,1,i+1,1)
!         res%d(i,1) = b%d(i,1) - mat%d(1,1)
!         res%d(i,1) = res%d(i,1)/this%d(i,i)
!       end do
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
      Integer(4) i, j
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
      q = da
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
    
    !В идеале нужно все циклы распараллелить.
    subroutine mtrx_hqr(this, q, r, maxr)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: q, r
      Integer(4), intent(in), optional :: maxr
      Type(Vector) c, piv, v, x
      Type(Mtrx) mat, mat2, da
      Integer(4) j, i, rank, info, Lwork
      Integer(4) n, m, maxrank
      DOUBLE PRECISION tmp, beta, tau
      DOUBLE PRECISION, Allocatable :: work(:), dtau(:)
      n = this%n
      m = this%m
      call c%init(m)
      call piv%init(m)
      do j = 1, m
        mat = this%submatrix(n, j, 1, j)
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
        
        tmp = piv%d(j)
        piv%d(j) = piv%d(rank)
        piv%d(rank) = tmp
        call da%swap(2, j, rank)
        tmp = c%d(j)
        c%d(j) = c%d(rank)
        c%d(rank) = tmp
        call x%init(n-rank+1)
        mat = da%submatrix(n,rank,rank,rank)
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
        mat = da%submatrix(n,m,rank,rank)
        !mat = mat - beta*mat2*((.T.mat2)*mat)
        call mat%update1(-1.0d0*beta,mat2,(.T.mat2)*mat)
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
      call r%pertcols(piv, -1)
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
    
    !Бесполезная. Или написать нормально, или через maxvol252.
    subroutine mtrx_geqr(this, q, r, maxr, pert)
      Class(Mtrx) :: this
      Type(Mtrx), intent(out) :: q, r
      Type(Vector), intent(inout), optional :: pert
      Integer(4), intent(in), optional :: maxr
      Type(Vector) c, w, piv, v, x
      Type(Mtrx) mat, mat2, da, dapr, AB, u1, u2, u
      Integer(4) j, i, rank, i1, j1, i2, j2
      Integer(4) n, m, maxrank, maxministeps, curministeps
      DOUBLE PRECISION tmp, beta, tau, ro, alpha, mu, f
      DOUBLE PRECISION, Allocatable :: dtau(:)
      
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
        mat = this%submatrix(n, j, 1, j)
        mat = (.T.mat) * mat
        c%d(j) = mat%d(1,1)
        piv%d(j) = j
      end do
      if (present(pert)) then
        piv = pert
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
        
        tmp = piv%d(j)
        piv%d(j) = piv%d(rank)
        piv%d(rank) = tmp
        call da%swap(2, j, rank)
        tmp = c%d(j)
        c%d(j) = c%d(rank)
        c%d(rank) = tmp
        call AB%swap(2, j, rank)
        call x%init(n-rank+1)
        mat = da%submatrix(n,rank,rank,rank)
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
        mat = da%submatrix(n,m,rank,rank)
        !mat = mat - beta*mat2*((.T.mat2)*mat)
        call mat%update1(-1.0d0*beta,mat2,(.T.mat2)*mat)
        do j = rank, m
          do i = rank, n
            da%d(i,j) = mat%d(i-rank+1,j-rank+1)
          end do
        end do
        do i = 2, n - rank + 1
          da%d(i+rank-1, rank) = 0.0d0
        end do
        mat = q%submatrix(n,n,1,rank)
        !mat = mat - (mat*(beta*mat2))*(.T.mat2)
        call mat%update1(-1.0d0*beta,mat*mat2,mat2)
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
          mat = AB%submatrix(rank-1,rank,1,rank)*da%submatrix(rank,m,rank,rank+1)/da%d(rank, rank)
        end if
        do j = 1, rank-1
          do i = rank+1, m
            AB%d(j, i) = AB%d(j,i) - mat%d(j, i-rank)
          end do
        end do
        if (rank > 1) then
          call mat%deinit()
        end if

        mat = AB%submatrix(rank,m,1,rank+1)
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
          tmp = piv%d(j1)
          piv%d(j1) = piv%d(j2)
          piv%d(j2) = tmp
          call da%swap(2, j1, j2)
          tmp = c%d(j1)
          c%d(j1) = c%d(j2)
          c%d(j2) = tmp
          call AB%swap(2, j1, j2)
          do i = i1+1, rank
            tmp = w%d(i-1)
            w%d(i-1) = w%d(i)
            w%d(i) = tmp
            tmp = piv%d(i-1)
            piv%d(i-1) = piv%d(i)
            piv%d(i) = tmp
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
            mat2 = mat * da%submatrix(i+1,m,i,i)
            do j = i, m
              da%d(i,j) = mat2%d(1,j-i+1)
              da%d(i+1,j) = mat2%d(2,j-i+1)
            end do
            mat2 = q%submatrix(n,i+1,1,i) * (.T.mat)
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
          mat = da%submatrix(n,rank,rank,rank)
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
          mat = da%submatrix(n,m,rank,rank)
          !mat = mat - beta*mat2*((.T.mat2)*mat)
          call mat%update1(-1.0d0*beta,mat2,(.T.mat2)*mat)
          do j = rank, m
            do i = rank, n
              da%d(i,j) = mat%d(i-rank+1,j-rank+1)
            end do
          end do
          do i = 2, n - rank + 1
            da%d(i+rank-1, rank) = 0.0d0
          end do
          mat = q%submatrix(n,n,1,rank)
          !mat = mat - (mat*(beta*mat2))*(.T.mat2)
          call mat%update1(-1.0d0*beta,mat*mat2,mat2)
          do j = rank, n
            do i = 1, n
              q%d(i,j) = mat%d(i,j-rank+1)
            end do
          end do
          call mat2%deinit()
          call mat%deinit()
          rank = rank - 1
          
          tmp = piv%d(rank)
          piv%d(rank) = piv%d(rank + 1)
          piv%d(rank + 1) = tmp
          dapr = 1.0d0*da
          call da%swap(2, rank, rank + 1)
          
          call givens(alpha, beta, da%d(rank,rank), da%d(rank+1,rank))
          call mat%init(2,2)
          mat%d(1,1) = alpha
          mat%d(1,2) = -beta
          mat%d(2,1) = beta
          mat%d(2,2) = alpha
          mat2 = mat * da%submatrix(rank+1,m,rank,rank)
          do j = rank, m
            da%d(rank,j) = mat2%d(1,j-rank+1)
            da%d(rank+1,j) = mat2%d(2,j-rank+1)
          end do
          mat2 = q%submatrix(n,rank+1,1,rank) * (.T.mat)
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
            mat = da%submatrix(rank-1,rank-1)
            u = mat%rtsolve(dapr%submatrix(rank-1,rank,1,rank))
          end if
          w%d(rank) = 1.0d0/da%d(rank,rank)/da%d(rank,rank)
          if (rank > 1) then
            u2 = AB%submatrix(rank-1,rank+1,1,rank+1) + AB%d(rank,rank+1)*u
          end if
          do i = 1, rank - 1
            w%d(i) = w%d(i) + u2%d(i,1)**2*w%d(rank) - u%d(i,1)**2/dapr%d(rank,rank)*2
          end do
          
          mu = AB%d(rank,rank+1)
          if (rank > 1) then
            u1 = AB%submatrix(rank-1,rank+1,1,rank+1)
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
            mat = u*(dapr%submatrix(rank,m,rank,rank+2)/dapr%d(rank,rank) - mu*da%submatrix(rank,m,rank,rank+2)/da%d(rank,rank))
            !mat = mat - u1*(da%submatrix(rank,m,rank,rank+2)/da%d(rank,rank))
            call mat%update1(-1.0d0/da%d(rank,rank),u1,da%submatrix(rank,m,rank,rank+2))
          end if
          do j = 1, rank-1
            do i = rank+2, m
              AB%d(j, i) = AB%d(j,i) + mat%d(j, i-rank-1)
            end do
          end do
          
          mat = AB%submatrix(rank,m,1,rank+1)
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
      call r%pertcols(piv, -1)
      q = q%submatrix(n, rank)
      call da%deinit()
      Deallocate(dtau)
      if (present(pert)) then
        pert = piv
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
      Integer(4), intent(in) :: n, m, r
      Integer(4) i, j
      Type(Mtrx) :: U, V
      DOUBLE PRECISION r1
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
      this%n = n
      this%m = m
      this%d = U%d
      call V%deinit()
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
  
    subroutine Mtrx_transform(this, vec)
      Class(Mtrx), intent(in) :: this
      Type(Vector), intent(out) :: vec
      Integer(4) i
      if (this%m == 1) then
        vec%n = this%n
        Allocate(vec%d(vec%n))
        do i = 1, vec%n
          vec%d(i) = this%d(i, 1)
        end do
      else
        vec%n = this%m
        Allocate(vec%d(vec%n))
        do i = 1, vec%n
          vec%d(i) = this%d(1, i)
        end do
      end if
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
    
    subroutine mtrx_replace(this, m2, coef, t, r, pert,pert2)
      Class(Mtrx) :: this
      Type(Mtrx) :: m2
      DOUBLE PRECISION, intent(in) :: coef
      Integer(4), intent(in) :: t
      Integer(4), intent(in) :: r
      Type(Vector), intent(in) :: pert, pert2
      DOUBLE PRECISION coef2
      Integer(4) i, j, m, tmp, tmp2
      if (t .eq. 1) then
        m = this%m
      else
        m = this%n
      end if
      if (coef .eq. 1.0d0) then
        do i = 1, r
          tmp = floor(pert%d(i) + 0.5)
          if (t == 1) then
            do j = 1, m
              tmp2 = floor(pert2%d(j) + 0.5)
              if (m2%d(tmp,tmp2) .ne. 0.0d0) then
                this%d(i,j) = m2%d(tmp,tmp2)
              end if
            end do
          else
            do j = 1, m
              tmp2 = floor(pert2%d(j) + 0.5)
              if (m2%d(tmp2,tmp) .ne. 0.0d0) then
                this%d(j,i) = m2%d(tmp2,tmp)
              end if
            end do
          end if
        end do
      else
        coef2 = coef - 1.0d0
        do i = 1, r
          tmp = floor(pert%d(i) + 0.5)
          if (t == 1) then
            do j = 1, m
              tmp2 = floor(pert2%d(j) + 0.5)
              if (m2%d(tmp,tmp2) .ne. 0.0d0) then
                this%d(i,j) = coef * m2%d(tmp,tmp2) - coef2 * this%d(i,j)
              end if
            end do
          else
            do j = 1, m
              tmp2 = floor(pert2%d(j) + 0.5)
              if (m2%d(tmp2,tmp) .ne. 0.0d0) then
                this%d(j,i) = coef * m2%d(tmp2,tmp) - coef2 * this%d(j,i)
              end if
            end do
          end if
        end do
      end if
    end
    
    function mtrx_tmnum(this, num) Result(res)
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
      !Совпадает 50/50. И возникает NaN.
      
      res%n = this%n
      res%m = this%m
      res%d = num * this%d
    end
    
    function mtrx_tmnumr(this, numin) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Real(4), intent(in) :: numin
      DOUBLE PRECISION num
      num = numin
      res%n = this%n
      res%m = this%m
      res%d = num * this%d
    end
    
    function mtrx_mnumt(num, this) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      DOUBLE PRECISION, intent(in) :: num
      res%n = this%n
      res%m = this%m
      res%d = num * this%d
    end
    
    function mtrx_mnumtr(num, this) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Real(4), intent(in) :: num
      res%n = this%n
      res%m = this%m
      res%d = num * this%d
    end
    
    function mtrx_divnum(this, num) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      DOUBLE PRECISION, intent(in) :: num
      res = this * (1.0d0 / num)
    end
    
    function mtrx_divnumr(this, num) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Real(4), intent(in) :: num
      res = this * (1.0d0 / num)
    end
    
    function mtrx_dotmul(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res
      res%n = this%n
      res%m = this%m
      if ((this%m == m2%m) .and. (this%n == m2%n)) then
        res%d = this%d * m2%d
      else
        print *, "error dotmul_mtrx"
      endif
    end
    
    function mtrx_dotvec(this, v) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Vector), intent(in) :: v
      Type(Mtrx) :: res
      Integer(4) i
      res%n = this%n
      res%m = this%m
      Allocate(res%d(res%n, res%m))
      if (this%n == v%n) then
        do i = 1, this%m
          res%d(:,i) = this%d(:,i) * v%d(:)
        end do
      else
        print *, "error dotvec_mtrx"
      endif
    end
    
    function mtrx_ddvec(this, v) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Vector), intent(in) :: v
      Type(Mtrx) :: res
      Integer(4) i
      res%n = this%n
      res%m = this%m
      Allocate(res%d(res%n, res%m))
      if (this%n == v%n) then
        do i = 1, this%m
          res%d(:,i) = this%d(:,i) / v%d(:)
        end do
      else
        print *, "error ddvec_mtrx"
      endif
    end
    
    function mtrx_dotdiv(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res
      res%n = this%n
      res%m = this%m
      if ((this%m == m2%m) .and. (this%n == m2%n)) then
        res%d = this%d / m2%d
      else
        print *, "error dotdiv_mtrx"
      endif
    end
  
    !Добавить отдельный оператор для транспонирования и проверять размеры, чтобы делать векторное умножение.
    function mtrx_mmtr(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx), intent(in) :: m2
      Type(Mtrx) :: res
      res%n = this%n
      res%m = m2%m
      if (this%m == m2%n) then
        !res%d = matmul(this%d, m2%d)
        Allocate(res%d(this%n, m2%m))
        call dgemm('N', 'N', this%n, m2%m, this%m, 1.0d0, this%d, this%n, m2%d, m2%n, 0.0d0, res%d, this%n)
      else
        print *, "error mul_mtrx", this%m, m2%n
      endif
    end
    
    function mtrx_tmvec(this, v) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Vector) :: res
      Type(Vector), intent(in) :: v
      res%n = this%n
      if (this%m == v%n) then
        !res%d = matmul(this%d, v%d)
        Allocate(res%d(this%n))
        call dgemv('N', this%n, this%m, 1.0d0, this%d, this%n, v%d, 1, 0.0d0, res%d, 1)
      else
        print *, "error mul_mvec"
      endif
    end
    
    function mtrx_mvect(v, this) Result(res)
      Type(Vector) :: res
      Type(Vector), intent(in) :: v
      Type(Mtrx), intent(in) :: this
      res%n = this%m
      if (this%n == v%n) then
        !res%d = matmul(v%d, this%d)
        Allocate(res%d(this%m))
        call dgemv('T', this%n, this%m, 1.0d0, this%d, this%n, v%d, 1, 0.0d0, res%d, 1)
      else
        print *, "error mul_mvec"
      endif
    end
    
    function mtrx_sum(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res
      if ((this%m == m2%m) .and. (this%n == m2%n)) then
        res%n = this%n
        res%m = this%m
        res%d = this%d + m2%d
      else
        print *, "error sum_mtrx"
      endif
    end
    
    function mtrx_transp(this) Result(res)
      Type(Mtrx), intent(in) :: this
      Type(Mtrx) :: res
      Integer(4) i, j
      res%n = this%m
      res%m = this%n
      Allocate(res%d(res%n, res%m))
      if (size(this%d,1) .ne. this%n) then
        print *, '1', size(this%d,1), this%n
      end if
      if (size(this%d,2) .ne. this%m) then
        print *, '2', size(this%d,2), this%m
      end if
      do i = 1, this%n
        do j = 1, this%m
          res%d(j,i) = this%d(i,j)
        end do
      end do
      !res%d = transpose(this%d(:,:this%m))
    end
    
    function mtrx_sub(this, m2) Result(res)
      Type(Mtrx), intent(in) :: this, m2
      Type(Mtrx) :: res
      if ((this%m == m2%m) .and. (this%n == m2%n)) then
        res%n = this%n
        res%m = this%m
        res%d = this%d - m2%d
      else
        print *, "error sub_mtrx"
        !call backtrace()
        stop
      endif
    end
    
    function tovec(this) Result(res)
      Type(Mtrx) :: this
      Type(vector) :: res
      call res%init(this%n*this%m)
      res%d = reshape(this%d, (/ this%n*this%m /))
    end
    
    function eye(n, m1) Result(res)
      Integer(4), intent(in) :: n
      Integer(4), intent(in), optional :: m1
      Type(Mtrx) :: res
      Integer(4) m, i
      if (present(m1)) then
        m = m1
      else
        m = n
      end if
      call res%init(n,m)
      do i = 1, min(m,n)
        res%d(i,i) = 1.0d0
      end do
    end
end
