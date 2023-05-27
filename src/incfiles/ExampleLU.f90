
subroutine ExampleLU()
  USE ModAppr
  Type(Mtrx) U, V, S, A, C, AR, E, Ahat, A1, AB, ABin, Q
  Type(IntVec) per1, per2, per3, peri
  Type(Vector) cin
  Integer(4) n, k, maxsteps, maxswaps, swapsmade, i, bigi, fulllu
  Double precision dsecnd, time
  Double precision avgf(7), avg2(7), avgc(7), avgt(7), maxf(7), max2(7), maxc(7), maxt(7), cur2
  
  !Desired size
  n = 1000
  !Desired rank
  k = 20
  !Maximum number of steps for maxvol and maxvol-rect
  !maxsteps = 1
  !Maximum number of row and column swaps for maxvol and maxvol-rect
  maxswaps = 100*4*k !Essentiall unlimited
  
  do bigi = 1, 1!00
  print *, 'Iteration', bigi
  
  write(*,'(A,I0,A,I0,A)') ' We create ', n ,' by ', n, ' matrix'
  print *, 'with random singular vectors (now generating)...'
  
!   RANDSVD ensemble
  call U%random(n)
  call V%random(n)
  call S%init(n,n)
  do i = 1, n
    S%d(i,i) = 1.0d0/2.0d0**i
  end do
  A = U*S*V
  
!   Block RANDSVD structure
!   call U%random(n/2,n/2)
!   call V%random(n/2,n/2)
!   call S%init(n/2,n/2)
!   do i = 1, n/2
!     S%d(i,i) = 1.0d0/2.0d0**i
!   end do
!   B1 = U*S*V
!   call U%random(n/2,n/2)
!   call V%random(n/2,n/2)
!   call S%init(n/2,n/2)
!   do i = 1, n/2
!     S%d(i,i) = 1.0d0/2.0d0**i
!   end do
!   B2 = U*S*V
!   call A%init(n,n)
!   A%d(1:n/2,1:n/2) = B1%d(1:n/2,1:n/2)
!   A%d(n/2+1:n,n/2+1:n) = B2%d(1:n/2,1:n/2)
!   call per1%permrand(n)
!   call per2%permrand(n)
!   call A%permrows(per1,1)
!   call A%permcols(per2,1)
!   call B1%gauss(n)
!   A = A + (B1*0.0000001d0)
!   call A%svd(U,S,V)
   
  call per1%perm(n)
  call per2%perm(n)
  
  print *, 'Done!'
  print *, ''
  
  write(*,'(A,I0,A)') ' We seek rank ', k, ' approximation, so Frobenius norm error of SVD is'
  cur2 = 0.0d0
  do i = k+1, n
    cur2 = cur2 + S%d(i,i)**2
  end do
  print *, 'SVD error:', sqrt(cur2)
  print *, 'SVD 2-norm error:', max(S%d(k+1,k+1), abs(S%d(k,k+1)))
  do i = 1, k
    S%d(i,i) = 0.0d0
    S%d(i+1,i) = 0.0d0
    S%d(i,i+1) = 0.0d0
  end do
  A1 = U*S*V
  print *, 'SVD C-norm error', A1%cnorm()
  print *, ''
  avgf(1) = avgf(1) + sqrt(cur2)
  avg2(1) = avg2(1) + max(S%d(k+1,k+1), abs(S%d(k,k+1)))
  avgc(1) = avgc(1) + A1%cnorm()
  maxf(1) = max(maxf(1), sqrt(cur2))
  max2(1) = max(max2(1), max(S%d(k+1,k+1), abs(S%d(k,k+1))))
  maxc(1) = max(maxc(1), A1%cnorm())
  
  !2) ADAPTIVE CROSS
  call per1%deinit()
  call per2%deinit()
  call per1%perm(n)
  call per2%perm(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  call ACA(n, n, k, A, C, AR)
  AR = .T.AR
  
  time = dsecnd() - time
  print *, '2) ACA time:', time

  E = A - C*AR
  print *, '2) ACA error:', E%fnorm()
  cur2 = E%norm2()
  print *, '2) ACA 2-norm error:', cur2
  print *, '2) ACA C-norm error:', E%cnorm()
  print *, ''
  
  avgf(2) = avgf(2) + E%fnorm()
  avg2(2) = avg2(2) + cur2
  avgc(2) = avgc(2) + E%cnorm()
  avgt(2) = avgt(2) + time
  maxf(2) = max(maxf(2), E%fnorm())
  max2(2) = max(max2(2), cur2)
  maxc(2) = max(maxc(2), E%cnorm())
  maxt(2) = max(maxt(2), time)
  
  !3) FULL GAUSS
  call per1%deinit()
  call per2%deinit()
  call per1%perm(n)
  call per2%perm(n)
  
  call E%copy(A)
  
  !Let's calculate the time
  time = dsecnd()
  
  call fgmaxvol(E, per1, per2, k)
  
  time = dsecnd() - time
  print *, '3) FULL GAUSS time:', time

  print *, '3) FULL GAUSS error:', E%fnorm()
  cur2 = E%norm2()
  print *, '3) FULL GAUSS 2-norm error:', cur2
  print *, '3) FULL GAUSS C-norm error:', E%cnorm()
  print *, ''
  
  avgf(3) = avgf(3) + E%fnorm()
  avg2(3) = avg2(3) + cur2
  avgc(3) = avgc(3) + E%cnorm()
  avgt(3) = avgt(3) + time
  maxf(3) = max(maxf(3), E%fnorm())
  max2(3) = max(max2(3), cur2)
  maxc(3) = max(maxc(3), E%cnorm())
  maxt(3) = max(maxt(3), time)
  
  !4) MAXVOL, 2 steps, r swaps
  maxsteps = 1
  call per1%deinit()
  call per2%deinit()
  call per1%perm(n)
  call per2%perm(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  !Find dominant k by k submatrix \hat A and construct C \hat A^{-1} R approximation
  call maxvol(Aelem, n, n, k, per1, per2, A, C, AR, maxsteps, k)
  
  time = dsecnd() - time
  print *, '4) MAXVOL LIMITED time:', time
  
  E = A - C*AR
  print *, '4) MAXVOL LIMITED error:', E%fnorm()
  cur2 = E%norm2()
  print *, '4) MAXVOL LIMITED 2-norm error:', cur2
  print *, '4) MAXVOL LIMITED C-norm error:', E%cnorm()
  print *, ''
  
  avgf(4) = avgf(4) + E%fnorm()
  avg2(4) = avg2(4) + cur2
  avgc(4) = avgc(4) + E%cnorm()
  avgt(4) = avgt(4) + time
  maxf(4) = max(maxf(4), E%fnorm())
  max2(4) = max(max2(4), cur2)
  maxc(4) = max(maxc(4), E%cnorm())
  maxt(4) = max(maxt(4), time)
  
  !5) RRLU by Pan
  !Pan, C.-T.: On the existence and computation of rank revealing LU factorizations. Linear Algebra Appl. 316, 199{222 (2000).
  !https://doi.org/10.1016/S0024-3795(00)00120-8
  call per1%deinit()
  call per1%perm(n)
  call per2%deinit()
  call per2%perm(n)

  call AR%deinit()
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
  !Let's limit ourselves to 100 swaps
  call A1%dominantr(2, k, n, per2, swapsmade, 100, Ahat, AB, ABin, cin)
  Ahat = Acols(Aelem, n, k, peri, per2, A)
  !Now we run maxvol in the chosen columns
  call Ahat%cmaxvol(per1, swapsmade, maxswaps, Q)
  
  call Q%permrows(per1, 2)
  AR = Arows(Aelem, k, n, per1, peri, A)

  time = dsecnd() - time
  print *, '5) RRLU time:', time
  
  E = A - Q*AR
  print *, '5) RRLU error:', E%fnorm()
  cur2 = E%norm2()
  print *, '5) RRLU 2-norm error:', cur2
  print *, '5) RRLU C-norm error:', E%cnorm()
  print *, ''
  
  avgf(5) = avgf(5) + E%fnorm()
  avg2(5) = avg2(5) + cur2
  avgc(5) = avgc(5) + E%cnorm()
  avgt(5) = avgt(5) + time
  maxf(5) = max(maxf(5), E%fnorm())
  max2(5) = max(max2(5), cur2)
  maxc(5) = max(maxc(5), E%cnorm())
  maxt(5) = max(maxt(5), time)
  
  !6) MAXVOL, unlimited steps
  maxsteps = 100
  call per1%deinit()
  call per2%deinit()
  call per1%perm(n)
  call per2%perm(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  !Find dominant k by k submatrix \hat A and construct C \hat A^{-1} R approximation
  call maxvol(Aelem, n, n, k, per1, per2, A, C, AR, maxsteps, maxswaps)
  
  time = dsecnd() - time
  print *, '6) UNLIMIT MAXVOL time:', time
  
  E = A - C*AR
  print *, '6) UNLIMIT MAXVOL error:', E%fnorm()
  cur2 = E%norm2()
  print *, '6) UNLIMIT MAXVOL 2-norm error:', cur2
  print *, '6) UNLIMIT MAXVOL C-norm error:', E%cnorm()
  print *, ''
  
  avgf(6) = avgf(6) + E%fnorm()
  avg2(6) = avg2(6) + cur2
  avgc(6) = avgc(6) + E%cnorm()
  avgt(6) = avgt(6) + time
  maxf(6) = max(maxf(6), E%fnorm())
  max2(6) = max(max2(6), cur2)
  maxc(6) = max(maxc(6), E%cnorm())
  maxt(6) = max(maxt(6), time)
  
  !9) OURS (3rho-locally maximum volume search), unlimited steps
  maxsteps = 100
  call per1%deinit()
  call per2%deinit()
  call per1%perm(n)
  call per2%perm(n)
  
  !Let's calculate the time
  time = dsecnd()
  
  call A1%copy(A)
  fulllu = 0
  !First run full Gaussian elimination
  call fgmaxvol(A1, per1, per2, k)
  !Since we have rho=1, which does not allow to limit the number of steps,
  !we can run maxvol first, since it at least does the swaps faster
  call maxvol(Aelem, n, n, k, per1, per2, A, C, AR, maxsteps, maxswaps)
  !Now we allow simultaneous exchanges of rows and columns
  call lumaxvol(Aelem, A, 1.0d0, per1, per2, k, fulllu, totsteps=swapsmade, cout=C, arout=AR)
  
  time = dsecnd() - time
  print *, '9) OURS time:', time
  !print *, swapsmade
  
  E = A - C*AR
  print *, '9) OURS error:', E%fnorm()
  cur2 = E%norm2()
  print *, '9) OURS 2-norm error:', cur2
  print *, '9) OURS C-norm error:', E%cnorm()
  print *, ''
  
  avgf(7) = avgf(7) + E%fnorm()
  avg2(7) = avg2(7) + cur2
  avgc(7) = avgc(7) + E%cnorm()
  avgt(7) = avgt(7) + time
  maxf(7) = max(maxf(7), E%fnorm())
  max2(7) = max(max2(7), cur2)
  maxc(7) = max(maxc(7), E%cnorm())
  maxt(7) = max(maxt(7), time)
  end do
  
  print *, avgf !Average Frobenius norm error
  print *, avg2 !Average spectral norm error
  print *, avgc !Average Chebyshev norm error
  print *, avgt !Average time
  print *, maxf !Maximum Frobenius norm error
  print *, max2 !Maximum spectral norm error
  print *, maxc !Maximum Chebyshev norm error
  print *, maxt !Maximum time
  
!   Output to file:
!   open(2, file = 'ExampleLresults')
!   write(2, *) avgf
!   write(2, *) avg2
!   write(2, *) avgc
!   write(2, *) avgt
!   write(2, *) maxf
!   write(2, *) max2
!   write(2, *) maxc
!   write(2, *) maxt
!   close(2)
end

!Gaussian elimination with complete pivoting
subroutine fgmaxvol(err, per1, per2, r, cura)
  USE ModMtrx
  Type(Mtrx) :: err
  Type(IntVec) :: per1, per2
  Integer(4) :: r
  Type(Mtrx), optional :: cura
  Type(Vector) :: u, v
  Integer(4) n, m, k, i, j, i1, ij1(2), ij2(2)
    
  n = err%n
  m = err%m
  k = r
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
    u = tovec(err%subarray(n, j, 1, j))
    v = tovec(err%subarray(i, m, i, 1))
    call err%update1v(-1.0d0/err%d(i,j),u,v)
    call per1%swap(i1, i)
    call per2%swap(i1, j)
    if (present(cura)) then
      call cura%swap(1, i1, i)
      call cura%swap(2, i1, j)
    end if
  end do
end

!rho-locally maximum volume search
!If full=1, also performs sorting to make the choice faster
!We only use full=0 (3rho-locally maximum volume search)
subroutine lumaxvol(Afun, Ain, ro, per1, per2, k, full, Ein, totsteps, cout, arout)
  USE ModAppr
      procedure(elem_fun) :: Afun !Function, returning elements of A
      Type(Mtrx), intent(in) :: Ain !Input matrix A (parameters for Afun)
      Double precision, intent(in) :: ro !Parameter rho
      Type(IntVec) :: per1, per2, persort !Permutations
      Integer(4), intent(in) :: k !Rank
      Integer(4), intent(in), optional :: full !If full=1, then rho-locally maximum volume is searched (instead of 3rho)
      Type(Mtrx), intent(in), optional :: Ein !Input error matrix, returned by Gaussian elimination
      Integer(4), intent(out), optional :: totsteps !Allowed number of steps
      Type(Mtrx), intent(out), optional :: cout, arout !Output low-rank factors $C$ and $\hat A^{-1} R$
      
      Type(Mtrx) A, C, R, CA, AR, Ahat, AI, AInew, E
      Type(Vector) CJ, CI, deltac, deltar, uvi1, uvj1, uvi2, uvj2, maxai, maxe
      Integer(4) n, m, i1, i2, j1, j2, i3, j3, ik, jk, i0, j0, curperm
      Integer(4) ij1(2), ij2(2)
      DOUBLE PRECISION CIJ, uc1, uc2
      Integer(4) fullused
      
      if (present(full)) then
        fullused = full
      else
        fullused = 0
      end if
      
      n = Ain%n
      m = Ain%m
      if (k > min(m,n)) then
        print *, "error in lumaxvol"
        print *, k, m, n
      end if
      A = Acols(Afun, n, m, per1, per2, Ain)
      C = A%subarray(n, k)
      R = A%subarray(k, n)
      Ahat = A%subarray(k, k)
      AI = .I.Ahat
      CA = C .dI. Ahat
      AR = Ahat .Id. R
      if (present(Ein)) then
        E = Ein
      else
        E = A - (C*AR)
      end if
      
      CIJ = ro
      
      !Check in CA
      ij1 = maxloc(CA%d)
      ij2 = minloc(CA%d)
      if (CA%d(ij1(1), ij1(2)) < -CA%d(ij2(1),ij2(2))) then
        ij1 = ij2
      end if
      if (abs(CA%d(ij1(1),ij1(2))) > CIJ) then
        i1 = ij1(2)
        i2 = ij1(1)
        j1 = i1
        j2 = i1
        CIJ = abs(CA%d(i2, i1) * AR%d(j1, j2) + AI%d(j1, i1) * E%d(i2, j2))
      end if
      
      !Check in AR
      ij1 = maxloc(AR%d)
      ij2 = minloc(AR%d)
      if (AR%d(ij1(1), ij1(2)) < -AR%d(ij2(1),ij2(2))) then
        ij1 = ij2
      end if
      if (abs(AR%d(ij1(1),ij1(2))) > CIJ) then
        j1 = ij1(1)
        j2 = ij1(2)
        i1 = j1
        i2 = j1
        CIJ = abs(CA%d(i2, i1) * AR%d(j1, j2) + AI%d(j1, i1) * E%d(i2, j2))
      end if
      
      !Fast check of E
      if ((fullused == 0) .and. (CIJ <= ro)) then
        ij1 = maxloc(E%d)
        ij2 = minloc(E%d)
        if (E%d(ij1(1), ij1(2)) > -E%d(ij2(1),ij2(2))) then
          i2 = ij1(1)
          j2 = ij1(2)
        else
          i2 = ij2(1)
          j2 = ij2(2)
        end if
        do i3 = 1, k
          do j3 = 1, k
            if (abs(CA%d(i2, i3) * AR%d(j3, j2) + AI%d(j3, i3) * E%d(i2, j2)) > CIJ) then
              i1 = i3
              j1 = j3
              CIJ = abs(CA%d(i2, i1) * AR%d(j1, j2) + AI%d(j1, i1) * E%d(i2, j2))
            end if
          end do
        end do
      end if
      
      !Full check
      if ((fullused == 1) .and. (CIJ <= ro)) then
        CJ = abs([CA%d])
        do i0 = 1, k
          CJ%d(i0+n*i0-n) = 0
        end do
        call persort%perm(n*k)
        call CJ%sort(persort)
        maxai = maxval(abs(AI%d),dim=2) !maxval(abs) is faster. But maxloc+minloc is faster than maxloc(abs)
        maxe = maxval(abs(E%d),dim=1)
        do j3 = k+1, m
          do i3 = 1, k
            do ik = n*k, k+1, -1
              jk = persort%d(ik)
              if (CJ%d(jk)*abs(AR%d(i3,j3)) + maxai%d(i3)*maxe%d(j3) <= CIJ) then
                exit
              end if
              i0 = 1 + (jk-1)/n
              j0 = 1 + mod(jk-1,n)
              if (abs(CA%d(j0,i0) * AR%d(i3,j3) + AI%d(i3,i0) * E%d(j0,j3)) > CIJ) then
                i1 = i0
                j1 = i3
                i2 = j0
                j2 = j3
                CIJ = abs(CA%d(i2, i1) * AR%d(j1, j2) + AI%d(j1, i1) * E%d(i2, j2))
              end if
            end do
          end do
        end do
        call persort%deinit()
        call CJ%deinit()
      end if
      
      curperm = 0
      do while (CIJ > ro)
        if (((i2 <= k) .and. (j2 <= k)) .or. (i1 > k) .or. (j1 > k)) then
          exit
        end if
        deltac = tovec(A%subarray(n, j2, 1, j2) - C%subarray(n, j1, 1, j1))
        deltar = tovec(A%subarray(i2, m, i2, 1) - R%subarray(i1, m, i1, 1))
        
        !Update AI from column of A
        CJ = tovec(AI%subarray(j1, k, j1, 1))
        call CI%init(k)
        CI%d = deltac%d(1:k)
        uc1 = -1.0d0/(1.0d0 + CI*CJ)
        uvi1 = AI*CI
        uvj1 = CJ
        call AI%update1v(uc1, uvi1, uvj1)
        call AInew%copy(AI)
        
        !Update AI from row of A
        CI = tovec(AI%subarray(k, i1, 1, i1))
        call CJ%init(k)
        CJ%d = deltar%d(1:k)
        CJ%d(j1) = A%d(i2,j2) - R%d(i1, j2)
        uc2 = -1.0d0/(1.0d0 + CI*CJ)
        uvi2 = CI
        uvj2 = CJ*AI
        call AI%update1v(uc2, uvi2, uvj2)
        
        !Update CA
        CI = tovec(AInew%subarray(j1, k, j1, 1))
        call E%update1v(-1.0d0, deltac, CI*R)
        call CA%update1v(1.0d0, deltac, CI)
        call E%update1v(-uc1, C*uvi1, uvj1*R)
        call CA%update1v(uc1, C*uvi1, uvj1)
        
        CJ = tovec(CA%subarray(n, i1, 1, i1))
        CI = tovec(CA%subarray(i2, k, i2, 1))
        CI%d(i1) = CI%d(i1) - 1.0d0
        call E%update1v(1.0d0/CA%d(i2, i1), CJ, CI*R)
        call CA%update1v(-1.0d0/CA%d(i2, i1), CJ, CI)
        
        CJ = tovec(CA%subarray(n, i1, 1, i1))
        call E%update1v(-1.0d0, CJ, deltar)
        
        !Update AR
        CJ = tovec(AR%subarray(j1, n, j1, 1))
        CI = tovec(AR%subarray(k, j2, 1, j2))
        CI%d(j1) = CI%d(j1) - 1.0d0
        call AR%update1v(-1.0d0/AR%d(j1, j2), CI, CJ)
        
        CI = tovec(AI%subarray(k, i1, 1, i1))
        call AR%update1v(1.0d0, CI, deltar)
        call AR%update1v(uc2, uvi2, uvj2*R)
        
        call A%swap(1, i1, i2)
        call A%swap(2, j1, j2)
        call E%swap(1, i1, i2)
        call E%swap(2, j1, j2)
        call CA%swap(1, i1, i2)
        call AR%swap(2, j1, j2)
        C = A%subarray(n, k)
        R = A%subarray(k, n)
        Ahat = A%subarray(k, k)
        
        call per1%swap(i1, i2)
        call per2%swap(j1, j2)
        
        CIJ = ro
        
        !Check in CA
        ij1 = maxloc(CA%d)
        ij2 = minloc(CA%d)
        if (CA%d(ij1(1), ij1(2)) < -CA%d(ij2(1),ij2(2))) then
          ij1 = ij2
        end if
        if (abs(CA%d(ij1(1),ij1(2))) > CIJ) then
          i1 = ij1(2)
          i2 = ij1(1)
          j1 = i1
          j2 = i1
          CIJ = abs(CA%d(i2, i1) * AR%d(j1, j2) + AI%d(j1, i1) * E%d(i2, j2))
        end if
        
        !Check in AR
        ij1 = maxloc(AR%d)
        ij2 = minloc(AR%d)
        if (AR%d(ij1(1), ij1(2)) < -AR%d(ij2(1),ij2(2))) then
          ij1 = ij2
        end if
        if (abs(AR%d(ij1(1),ij1(2))) > CIJ) then
          j1 = ij1(1)
          j2 = ij1(2)
          i1 = j1
          i2 = j1
          CIJ = abs(CA%d(i2, i1) * AR%d(j1, j2) + AI%d(j1, i1) * E%d(i2, j2))
        end if
        
        !Fast check of E
        if ((fullused == 0) .and. (CIJ <= ro)) then
          ij1 = maxloc(E%d)
          ij2 = minloc(E%d)
          if (E%d(ij1(1), ij1(2)) > -E%d(ij2(1),ij2(2))) then
            i2 = ij1(1)
            j2 = ij1(2)
          else
            i2 = ij2(1)
            j2 = ij2(2)
          end if
          CIJ = 0
          do i3 = 1, k
            do j3 = 1, k
              if (abs(CA%d(i2, i3) * AR%d(j3, j2) + AI%d(j3, i3) * E%d(i2, j2)) > CIJ) then
                i1 = i3
                j1 = j3
                CIJ = abs(CA%d(i2, i1) * AR%d(j1, j2) + AI%d(j1, i1) * E%d(i2, j2))
              end if
            end do
          end do
        end if
        
        !Full check
        if ((fullused == 1) .and. (CIJ <= ro)) then
          CJ = abs([CA%d])
          do i0 = 1, k
            CJ%d(i0+n*i0-n) = 0
          end do
          call persort%perm(n*k)
          call CJ%sort(persort)
          maxai = maxval(abs(AI%d),dim=2)
          maxe = maxval(abs(E%d),dim=1)
          do j3 = k+1, m
            do i3 = 1, k
              do ik = n*k, k+1, -1
                jk = persort%d(ik)
                if (CJ%d(jk)*abs(AR%d(i3,j3)) + maxai%d(i3)*maxe%d(j3) <= CIJ) then
                  exit
                end if
                i0 = 1 + (jk-1)/n
                j0 = 1 + mod(jk-1,n)
                if (abs(CA%d(j0,i0) * AR%d(i3,j3) + AI%d(i3,i0) * E%d(j0,j3)) > CIJ) then
                  i1 = i0
                  j1 = i3
                  i2 = j0
                  j2 = j3
                  CIJ = abs(CA%d(i2, i1) * AR%d(j1, j2) + AI%d(j1, i1) * E%d(i2, j2))
                end if
              end do
            end do
          end do
          call persort%deinit()
          call CJ%deinit()
        end if
      
        curperm = curperm + 1
      end do
      
      if (present(totsteps)) then
        totsteps = curperm
      end if
      if (present(cout)) then
        cout = C
        arout = AR
        call cout%permrows(per1, 2)
        call arout%permcols(per2, 2)
      end if
    end
    

!Adaptive cross approximation aka Cross 2D
!Written to stop at rank MaxRank instead of based on the error bound
!Initial code provided by Stanislav Stavtsev, see, e.g.:
!Aparinov, A.A., Setukha, A.V. & Stavtsev, S.L.
!Parallel Implementation for Some Applications of Integral Equations Method. 
!Lobachevskii J Math 39, 477–485 (2018). https://doi.org/10.1134/S1995080218040029
subroutine ACA(Ni, Nj, MaxRank, param, U, V)
     USE ModAppr
       Integer, intent(in) :: Ni, Nj, MaxRank !Matrix sizes, desired rank
       Type(Mtrx), intent(in) :: param !Matrix parameters
       Type(Mtrx), intent(out) :: U, V !Output low-rank factors
       
       Logical, allocatable :: KnR(:), KnC(:) !Used rows and columns
       Integer(4) :: iP, jP !Old row and column indices
       Integer(4) :: jN !New column index
       Integer(4) :: NSkel !Current skeletion (submatrix) size
       Double precision :: Elem, ElemS !Pivot element and its square root
       Double precision :: MaxElC, MaxElR !Max elements in column and row
       Integer(4) i, j, s

        call U%init(Ni,maxrank)
        call V%init(Nj,maxrank)
        Allocate(KnR(Ni))
        Allocate(KnC(Nj))

        KnR = .false.
        KnC = .false.
        
        !Init column indices 
        jN = 1
        jP = 1
        
        NSkel  = 0

        !Start approximation process
        do while (NSkel < MaxRank)
            !Number of skeletons to aproximate matrix
            NSkel = NSkel + 1                    

            !Calculate elements in column
            j  = jP
            do i = 1, Ni
              Elem = Aelem(i,j,param)             
              if (.not. KnR(i)) then             
                do s = 1, NSkel - 1
                  Elem = Elem - U%d(i,s) * V%d(j,s)
                end do
                U%d(i,NSkel) = Elem
              else
                U%d(i,NSkel) = 0
              end if
            end do

            !Calculate the maximum element in column
            MaxElC = 0
            iP = 0
            do i = 1, Ni
              if ((.not. KnR (i)) .and. (abs(U%d(i,NSkel)) >= MaxElC)) then
                MaxElC = abs(U%d(i,NSkel))
                iP = i
              end if
            end Do
            
            !Calculate elements in row
            i = iP
            do j = 1, Nj
              if (.not. KnC(j)) then 
                Elem = Aelem(i,j,param)
                do s = 1, NSkel - 1
                  Elem = Elem - U%d(i,s) * V%d(j,s)
                end do
                V%d(j,NSkel) = Elem
              else
                V%d(j,NSkel) = 0
              end if
            end do

            !Calculate the maximum element in row
            MaxElR = 0
            jN = 0
            do j = 1, Nj
              if ((.not. KnC(j)) .and. (abs(V%d(j,NSkel)) >= MaxElR) .and. (j .ne. jP)) then
                MaxElR = abs(V%d(j,NSkel))
                jN = j
              End If
            End Do
!      ..  
!      ..  Form skeletons at cross iP, jP
!      ..  U part of skeleton is already ready. All we shall do is to
!      ..  norm vector V(jP,NSkel) = MaxElC to 1. This value is Elem. If Elem
!      ..  is equal to zero, i.e. MaxElC = 0, then U is changed, not V
!      ..  Element U(iP,NSkel) is changed to 1 in this case.
!      ..
            
            Elem = U%d(iP,NSkel)
            ElemS = sqrt(abs(Elem))
            if (MaxElC > 0) then
                do j = 1, Ni
                   U%d(j,NSkel) = ElemS*U%d(j,NSkel)/Elem
                end do
                do j = 1, Nj
                   V%d(j,NSkel) = V%d(j,NSkel)/ElemS
                end do
            else
                U%d(iP,NSkel) = 1.0d0
            end if 

            !Mark new row and column
            KnR (iP) = .true.
            KnC (jP) = .true.
            jP = jN

        end do

        Deallocate(KnR,KnC)
        
     end subroutine
