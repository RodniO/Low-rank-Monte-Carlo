Module ModRecon
  USE ModSparse
  
  !Module for matrix reconstruction
  
  Contains

  !Reconstruction subroutine. The only one to be called externally. Others just split the code.
  subroutine reconstruct(KnownSparse, r, lcurmat, rcurmat, maxtime, maxerror, verbose, maxio)
    Type(SparseRow), intent(in) :: KnownSparse !Known elements as a sparse matrix
    Integer(4), intent(in) :: r !Desired rank
    Type(Mtrx) :: lcurmat, rcurmat !Low-rank factors. Must be initialized in advance
    Double precision, intent(in) :: maxtime !Maximum allowed computation time
    Double precision, optional, intent(in) :: maxerror !Maximum allowed relative error
    Integer(4), optional, intent(in) :: verbose !How much to print. 0-3.
    !0 - no print
    !1 - no explanation
    !2 - with explanation
    !3 - no explanation, with convergence info
    Integer(4), optional, intent(out) :: maxio
  
    Double precision dsecnd !Function for time evaluation from lapack
    Double precision maxerror_ !Absolute allowed maximum error
    Double precision curerror !Temporary for error computation
    Double precision curtime !Current time
    Integer(4) n, m !Number of rows and columns
    Integer(4) curr !Current rank
  
    Type(SparseCol) KnownSparsec !CSC version of KnownSparse
    Type(IntVec) perm1, perm2 !Permutations of rows and columns
    Type(Mtrx) mask !Partial mask of nonzero elements
    Type(SparseRow) smask !Sparse mask (zeros and ones)
    Integer(4) r1 !Initial approximation rank
    Integer(4) r2 !Maximum number of rows and columns
    Double precision q !Sparsity fraction
    Double precision sq !Step size
    Double precision rtransl !mu-coherence bound. Matrices must have mu << rtransl^2
    Integer(4) prevbads, maxbads !Current and maximum number of bad steps per rank increase
    Double precision matmasknorm !Total norm of known elements
    Double precision totmult !Convergence multiplier estimate
    Integer(4) maxi !Total number of iterations
    Integer(4) maxswaps !Maximum number of maxvol swaps
    Type(Mtrx) p1, p2, cpmat, rpmat, u, s, v !Temporary matrices
    Double precision tmp !Temporary
    Integer(4) i, j, j1, j2, k, l, nz !Cycle counters
    Integer(4) verbose_ !Another name for verbose
  
    !Variables for iteration cycle (see below)
    Integer(4) everybi, bi, bads, goods, steps, swaps, maxsteps, lastr, maxir
    Double precision minerror, multmin, dsq, minsq, cerror, perror
  
    !Names for saved variables
    Integer(4) savedi, savedbi, savedr, savedmb
    Type(Mtrx) savedlcm, savedrcm, savedcpm, savedrpm
    Type(IntVec) savedp1, savedp2
    
    verbose_ = 0
    if (present(verbose)) then
      verbose_ = verbose
    end if
  
    KnownSparsec = KnownSparse
  
    n = lcurmat%n
    m = rcurmat%m
    curr = lcurmat%m
    q = dble(KnownSparse%nz)/n/m
  
    !Initializing parameters
    !CAN BE CHANGED IF NECESSARY
    prevbads = 5 !Larger value means longer rank increase
    maxbads = prevbads
    r1 = 2*r
    r2 = r1 + floor(0.7d0*r/q + 0.99d0) !0.7 is from experiments. Increase/decrease if needed.
    sq = 0.5d0 / q !0.75/q is theoretically optimal
    rtransl = 2.0d0*sqrt(log(2.0d0*n*r)) !Lower value means longer rank increase
  
    r2 = min(r2, n)
    r2 = min(r2, m)
    sq = max(sq, 1.0d0)
    
    if (verbose_ > 0) then
      if (verbose_ == 2) then
        print *, 'Size; rank; maxvol rank; rows and columns in total.'
      end if
      write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') ' ', n, ' x ', m, ', ', r, ', ', r1, ', ', r2
    end if
    
    curtime = dsecnd()
  
    matmasknorm = KnownSparse%fnorm()
    if (present(maxerror)) then
      maxerror_ = maxerror*matmasknorm
    else
      maxerror_ = eps*matmasknorm
    end if
  
    call perm1%perm(n)
    call perm2%perm(m)
  
    !Search for submatrix with nonzero diagonal (thus usually full rank)
    mask = KnownSparse%detransform(r1)
    do i = 1, r1
      if (mask%d(i,i) == 0.0d0) then
        do l = i, r1
          do k = i, m
            if (mask%d(l,k) .ne. 0.0d0) then
              call mask%swap(2, i, k)
              call mask%swap(1, i, l)
              call perm1%swap(i, l)
              call perm2%swap(i, k)
              exit
            end if
          end do
          if (mask%d(i,i) .ne. 0.0d0) then
            exit
          end if
        end do
      end if
    end do

    !Initializing saved variables
    savedi = 0
    savedbi = 0
    savedr = curr
    call savedlcm%copy(lcurmat)
    call savedrcm%copy(rcurmat)
    call savedcpm%copy(lcurmat)
    call savedrpm%copy(rcurmat)
    call savedp1%copy(perm1)
    call savedp2%copy(perm2)
    savedmb = maxbads
    
    !Initializing cycle variables
    everybi = ((n+m)*r)/r2**2/2+1 !How often to check the error
    minerror = matmasknorm !Minimum error reached
    multmin = 0.0d0 !Minimum convergence multiplier
    minsq = sq !Minimum step size reached
    dsq = 0.0d0 !Step size increase
    bi = 0 !Number of approximations
    i = 0 !Number of steps
    lastr = 0 !Maximum rank reached
    bads = 0 !'Bad' steps
    goods = 0 !'Good' steps
    if (present(maxio)) then
      maxio = 0 !Output of maximum steps
    end if
    steps = 0 !Maxvol steps
    swaps = 0 !Maxvol swaps
    maxsteps = 10 !Maximum number of maxvol steps
    maxir = 100000000 !Steps until full rank. Initialized as infinity
    perror = matmasknorm**2 !Previous step error in the final SVD
    cerror = 0.5d0*perror !Current step error in the final SVD
    do while (1 > 0)
      i = i + 1
    
      !Evaluate maximum i, after which to stop rank increase
      !and check for divergence
      if ((i/10)*10 == i) then
        if (dsecnd() > curtime) then
          maxir = floor(0.5d0/((dsecnd()-curtime)/maxtime/i))
        end if
        if (rcurmat%fnorm() > 100.0d0*matmasknorm/sqrt(q)) then
          print *, 'Divergence occured!'
          smask = KnownSparse%mask()
          p1 = ((lcurmat*rcurmat) .dot. smask%detransform()) - KnownSparse%detransform()
          print *, 'Reached relative error:', p1%fnorm()/matmasknorm
          exit
        end if
      end if
      
      !Increase sizes of the low-rank factors
      !if rank was increased
      if (curr > lcurmat%m) then
        call p2%init(n,1)
        lcurmat = lcurmat%unite(2,p2)
        call p2%deinit()
        call p2%init(1,m)
        rcurmat = rcurmat%unite(1,p2)
        call p2%deinit()
      end if
      call cpmat%copy(lcurmat)
      call rpmat%copy(rcurmat)
      cpmat%m = curr
      rpmat%n = curr
      if (curr > lcurmat%m) then
        call lcurmat%copy(cpmat)
        call rcurmat%copy(rpmat)
      end if
      
      !Choose maximum number of steps for maxvol
      if (maxsteps > 0) then
        maxswaps = r
      end if
      if (i*2 >= maxir) then
        maxsteps = 2
      else
        maxsteps = 10
      end if
      
      !Check if we can stay in the same rows and columns
      totmult = cerror/perror
      if (maxswaps > 0) then
        totmult = 1.0d0
      end if
      multmin = min(multmin, totmult)
      !if (1.0d0 - totmult < 0.5d0*(1.0d0-multmin)) then
      if (1.0d0 - totmult < 2.0d0*(1.0d0-multmin)) then
        maxswaps = r
        multmin = 1.0d0
      else
        maxswaps = 0
      end if
      if (maxswaps == 0) then
        maxsteps = 0
      end if
      
      !Run approximation
      perror = cerror
      call ApproxSparse(lcurmat, rcurmat, KnownSparse, KnownSparsec, &
                curr, r1, r2, perm1, perm2, &
                steps, swaps, maxsteps, maxswaps, sq, cerror)
      
      !Step choice.
      !Saves data and reloads it in case of divergence
      if (curr <= r) then
        if (maxsteps .ne. 0) then
          bi = bi + 1
          if (((bi/everybi)*everybi == bi) .or. (rcurmat%fnorm() > 10.0d0*matmasknorm/sqrt(q))) then
            if (maxsteps .ne. 0) then
              curerror = 0.0d0
              j = 0
              do nz = 1, KnownSparse%nz
                do while (KnownSparse%i(j+1) <= nz)
                  j = j + 1
                end do
                j2 = KnownSparse%j(nz)
                tmp = -KnownSparse%d(nz)
                do j1 = 1, lcurmat%m
                  tmp = tmp + lcurmat%d(j,j1)*rcurmat%d(j1,j2)
                end do
                curerror = curerror + tmp**2
              end do
              curerror = sqrt(curerror)
              if (verbose_ == 3) then
                print *, curerror
              end if
            end if
            if (curerror > 2.0d0*minerror) then
              sq = (sq - 1.0d0)/2.0d0
              sq = max(sq, 1.0d0)
              dsq = sq/2.0d0
              sq = min(sq,0.5d0/q)
              minsq = sq
              dsq = min(dsq, (0.5d0/q-sq)/2.0d0)
              if (present(maxio)) then
                maxio = maxio + i - savedi
              end if
              i = savedi
              bi = savedbi
              curr = savedr
              call lcurmat%copy(savedlcm)
              call rcurmat%copy(savedrcm)
              call cpmat%copy(savedcpm)
              call rpmat%copy(savedrpm)
              call perm1%copy(savedp1)
              call perm2%copy(savedp2)
              maxbads = savedmb
              cerror = perror
              multmin = 0.0d0
              prevbads = maxbads
              bads = 0
              goods = 0
            else
              if (curerror < minerror) then
                minerror = curerror
                sq = sq + dsq
                dsq = dsq/2.0d0
                savedi = i
                savedbi = bi
                savedr = curr
                call savedlcm%copy(lcurmat)
                call savedrcm%copy(rcurmat)
                call savedcpm%copy(cpmat)
                call savedrpm%copy(rpmat)
                call savedp1%copy(perm1)
                call savedp2%copy(perm2)
                savedmb = maxbads
              end if
            end if
          end if
        end if
      end if
      
      !Perform rank increase if necessary
      if (curr < r) then
        cpmat = cpmat%unite(2, (-1.0d0)*lcurmat)
        rpmat = rpmat%unite(1, (-1.0d0)*rcurmat)
        call cpmat%qr(u, s)
        call rpmat%lq(p2,v)
        s = s * p2
        call s%svd(cpmat, p2, rpmat)
        call p2%deinit()
        u = u*cpmat
        if ((i*r > curr*maxir) .and. (curr < r)) then
          curr = curr + 1
          cerror = perror
          multmin = 0.0d0
          prevbads = maxbads
          bads = 0
          goods = 0
        else if (u%cnorm()*sqrt(1.0d0*n) < rtransl) then
          goods = goods + 1
          if (goods > prevbads) then
            curr = curr + 1
            cerror = perror
            multmin = 0.0d0
            prevbads = maxbads
            bads = 0
            goods = 0
          end if
        else if (u%cnorm()*sqrt(1.0d0*n) < rtransl) then
          bads = 0
        else if (u%cnorm()*sqrt(1.0d0*n) > rtransl) then
          bads = bads + 1
          maxbads = max(bads, maxbads)
        end if
      end if
 
      !Try to use larger step if rank increase is over
      if ((i > maxir) .and. (lastr == 0)) then
        dsq = (0.5/q - sq)/2.0d0
        lastr = 1
      end if
      
      !Exit if time is up or desired quality reached
      if ((dsecnd() >= maxtime + curtime) .or. (minerror <= maxerror_)) then
        maxi = i
        exit
      end if

    end do
    
    if (present(maxio)) then
      maxio = maxio + maxi
    end if

    !Debug output
    if (verbose_ == 2) then
      print *, 'Iterations, maxvol steps, maxvol swaps'
    end if
    if (verbose_ > 0) then
      print *, maxi, steps, swaps
    end if
    if (verbose_ == 3) then
      print *, (steps+0.0d0)/i, (swaps+0.0d0)/maxi/2, sq, minsq
    end if
    
  end
  
  !Replaces elements of the matrix according to sparse mask and step coef.
  !For internal calls ONLY.
  subroutine replace(this, srmat, scmat, coef, t, r, perm, perm2)
    Type(Mtrx) :: this
    Type(SparseRow), intent(in) :: srmat
    Type(SparseCol), intent(in) :: scmat
    DOUBLE PRECISION, intent(in) :: coef
    Integer(4), intent(in) :: t
    Integer(4), intent(in) :: r
    Type(IntVec), intent(in) :: perm, perm2
    DOUBLE PRECISION coef2
    Integer(4) i, j, nz, m, k
    Type(IntVec) perm2i

    if (t .eq. 1) then
      m = this%m
    else
      m = this%n
    end if
    perm2i = perm2%perminv()
    if (coef .ne. 1) then
      coef2 = coef - 1.0d0
      if (t == 1) then
        do i = 1, r
          j = perm%d(i)
          do nz = srmat%i(j), srmat%i(j+1)-1
            k = perm2i%d(srmat%j(nz))
            this%d(i,k) = coef * srmat%d(nz) - coef2 * this%d(i,k)
          end do
        end do
      else
        do i = 1, r
          j = perm%d(i)
          do nz = scmat%j(j), scmat%j(j+1)-1
            k = perm2i%d(scmat%i(nz))
            this%d(k,i) = coef * scmat%d(nz) - coef2 * this%d(k,i)
          end do
        end do
      end if
    else
      if (t == 1) then
        do i = 1, r
          j = perm%d(i)
          do nz = srmat%i(j), srmat%i(j+1)-1
            k = perm2i%d(srmat%j(nz))
            this%d(i,k) = srmat%d(nz)
          end do
        end do
      else
        do i = 1, r
          j = perm%d(i)
          do nz = scmat%j(j), scmat%j(j+1)-1
            k = perm2i%d(scmat%i(nz))
            this%d(k,i) = scmat%d(nz)
          end do
        end do
      end if
    end if
  end

  !Approximation subroutine for reconstruction.
  !For internal calls ONLY.
  subroutine ApproxSparse(lcurmat, rcurmat, KnownSparse, KnownSparsec, &
                r, r1, r2, perm1, perm2, &
                steps, swaps, maxsteps, maxswaps, sq, error)
    Type(Mtrx) :: lcurmat, rcurmat
    Type(SparseRow), intent(in) :: KnownSparse
    Type(SparseCol), intent(in) :: KnownSparsec
    Type(Mtrx) u, s, v, cs, us, rs, cs1, rs1, curc
    Type(Mtrx) curc2, curr2
    Type(IntVec) :: perm1, perm2
    Type(Vector) tau1, tau2
    Type(Vector) work
    Integer(4), intent(in) :: r, r1, r2
    Integer(4) :: steps, swaps, maxsteps
    Integer(4) :: maxswaps
    Double precision, intent(in) :: sq
    Double precision, intent(out) :: error
    Integer(4) st, tst, n, prevst, prevtst
    Integer(4) i, j, j1
    
    n = lcurmat%n
    call lcurmat%permrows(perm1, 1)
    call rcurmat%permcols(perm2, 1)
    
    st = 0
    prevst = 1
    prevtst = swaps
    do i = 1, maxsteps/2
      curc = lcurmat * rcurmat%subarray(r,r1)
      call replace(curc, KnownSparse,KnownSparsec,sq,2,r1,perm2,perm1)
      call curc%cmaxvol(perm1, tst, r1, curc2, lcurmat, 1)
      swaps = swaps + tst
      st = st + 1
      if ((prevst .eq. 0) .and. (tst .eq. 0)) exit
      prevst = tst
      curc = lcurmat%subarray(r1,r) * rcurmat
      call replace(curc, KnownSparse,KnownSparsec,sq,1,r1,perm1,perm2)
      curc = .T.curc
      call curc%cmaxvol(perm2, tst, r1, curr2, rcurmat, 2)
      swaps = swaps + tst
      st = st + 1
      if ((prevst .eq. 0) .and. (tst .eq. 0)) exit
      prevst = tst
    end do
    steps = steps + st
    
    if (maxswaps > 0) then
      if (tst > 0) then
        curc = lcurmat * rcurmat%subarray(r,r1)
        call replace(curc,KnownSparse,KnownSparsec,sq,2,r1,perm2,perm1)
        curc2 = curc .dI. curc%subarray(r1, r1)
      end if
      call curc2%hmaxvol2(1, r1, r2, perm1, .true., lcurmat)
    end if

    if (maxswaps > 0) then
      curr2 = .T.curr2
      call curr2%hmaxvol2(2, r1, r2, perm2, .true., rcurmat)
    end if
    
    if (swaps - prevtst >= r) then
      maxswaps = r
    else
      maxswaps = 0
    end if
    
    cs = lcurmat * rcurmat%subarray(r,r2)
    call replace(cs, KnownSparse,KnownSparsec,sq,2,r2,perm2,perm1)
    us = cs%subarray(r2, r2)
    call us%svd(u, s, v)
    u = u%subarray(r2, r1)
    s = s%subarray(r1, r1)
    v = v%subarray(r1, r2)
    do j = 1, r1
      s%d(j,j) = 1.0d0/s%d(j,j)
    end do
    cs = cs * (.T.v)
    us = s*(.T.u)
    rs = lcurmat%subarray(r2,r) * rcurmat
    call replace(rs, KnownSparse,KnownSparsec,sq,1,r2,perm1,perm2)
    rs = us*rs
    
    call cs%halfqr(u, tau1, cs1)
    call rs%halflq(rs1, tau2, v)
    
    us = cs1 * rs1
    call us%svd(cs, s, rs)
    cs = cs%subarray(r1, r)
    error = 0.0d0
    do i = r+1, r1
      error = error + s%d(i,i)**2
    end do
    s = s%subarray(r, r)
    rs = rs%subarray(r, r1)

    rs = s*rs
    
    cs = cs%multq(u, tau1, 'L', 'D')
    rs = rs%multq(v, tau2, 'R', 'U')
    
    call cs%permrows(perm1, 2)
    call rs%permcols(perm2, 2)

    if (maxsteps < 0) then
      call u%deinit()
      call us%deinit()
      call v%deinit()
      call u%init(KnownSparse%nz, r)
      do i = 1, r
        j1 = 0
        do j = 1, u%n
          do while (KnownSparse%i(j1+1) <= j)
            j1 = j1 + 1
          end do
          u%d(j,i) = cs%d(j1,i) * rs%d(i,KnownSparse%j(j))
        end do
      end do
      call v%init(u%n, 1)
      do i = 1, u%n
        v%d(i,1) = KnownSparse%d(i)
      end do
      call work%init(2*u%n*u%m**2)
      call dgels('N', u%n, u%m, 1, u%d, u%n, v%d, u%n, work%d, 2*u%n*u%m**2, j1)
      do i = 1, r
        s%d(i,i) = v%d(i,1)
      end do   
    end if 
    
    lcurmat = cs
    rcurmat = rs
  end

end module
