
!Monte-Carlo examples.
subroutine ExampleM()
  USE ModMtrx
  USE ModMonte
  
  Type(Vector) temps
  Type(IntVec) nv
  Integer(4) i
  Double precision n0
  
  Integer(4) n
  Double precision np, tmpi
  
  !call init_random_seed()
  
  n = 10
  
  print *, ''
  print *, 'CLASSICAL SMOLUCHOWSKI EQUATIONS'
  
  call nv%init(n)
  nv%d(1) = 100000
  n0 = 1.0d0
  print *, ''
  print *, 'Wagner (majorant) method:'
  call MonteWagnerSimple(n0, nv, 3.0d0, 1)
  call nv%deinit()
  call nv%init(n)
  nv%d(1) = 100000
  n0 = 1.0d0
  print *, ''
  print *, 'Walker-Sabelfeld (stratified) method:'
  call MonteWalker(n0, nv, 3.0d0, 1)
  call nv%deinit()
  call nv%init(n)
  nv%d(1) = 100000
  n0 = 1.0d0
  print *, ''
  print *, 'Low-rank method:'
  call MonteSimple(n0, nv, 3.0d0, 1)
  call nv%deinit()
  call nv%init(n)
  nv%d(1) = 100000
  n0 = 1.0d0
  print *, ''
  print *, 'Low-rank method with array of particles:'
  call MonteSimplePart(n0, nv, 3.0d0, 1) !Use when N_p (number of particles) << M (maximum particle mass)
  call nv%deinit()
  call nv%init(n)
  nv%d(1) = 100000
  n0 = 1.0d0
  print *, ''
  print *, 'Acceptance-rejection method:'
  call MonteUpB(n0, nv%d(1)/10, 3.0d0, 1)
  call nv%deinit()
  call nv%init(n)
  nv%d(1) = 100000
  n0 = 1.0d0
  print *, ''
  print *, 'Inverse method:'
  call MonteInverse(n0, nv%d(1)/10, 3.0d0, 1)
  call nv%deinit()
  
  print *, ''
  print *, 'TEMPERATURE-DEPENDENT EQUATIONS'
  
  call nv%init(n)
  call temps%init(n)
  nv%d(1) = 10000
  temps%d = 1
  n0 = 0.3d0/pi
  print *, ''
  print *, 'Low-rank method:'
  call MonteTemp(n0, nv, temps, 1000.0d0, 0.4d0, 0.1d0, 1)
!   call nv%deinit()
!   call temps%deinit()
!   call nv%init(n)
!   call temps%init(n)
!   nv%d(1) = 10000
!   temps%d = 1
!   n0 = 0.3d0/pi
!   print *, ''
!   print *, 'Wagner (majorant) method:'
!   call MonteWagnerRanked(3, n0, nv, temps, 1000.0d0, 0.4d0, 0.1d0, 1) !INCORRECT
  call nv%deinit()
  call temps%deinit()
  call nv%init(n)
  call temps%init(n)
  nv%d(1) = 10000
  temps%d = 1
  n0 = 0.3d0/pi
  print *, ''
  print *, 'Acceptance-rejection method:'
  call MonteUpBtemp(n0, nv%d(1), temps%d(1), 1000.0d0, 1)
  call nv%deinit()
  call temps%deinit()
  
  print *, ''
  print *, 'BOLTZMANN EQUATIONS'
  
  call nv%init(2**10)
  call temps%init(2**10)
  nv%d(1) = 999
  i = 1
  tmpi = 1.0d0
  np = nv%d(1)
  do while (tmpi > 0)
    i = i + 1
    tmpi = floor(nv%d(1)*exp(-0.01d0*(i-1)))
    nv%d(i) = floor(nv%d(1)*exp(-0.01d0*(i-1)))
    np = np + tmpi
  end do
  temps%d(:) = 1
  n0 = 0.1d0
  print *, ''
  print *, 'Low-rank method:'
  call MonteDSMC(n0, nv, temps, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1)
  call nv%deinit()
  call temps%deinit()
  call nv%init(2**10)
  call temps%init(2**10)
  nv%d(1) = 999
  i = 1
  tmpi = 1.0d0
  np = nv%d(1)
  do while (tmpi > 0)
    i = i + 1
    tmpi = floor(nv%d(1)*exp(-0.01d0*(i-1)))
    nv%d(i) = floor(nv%d(1)*exp(-0.01d0*(i-1)))
    np = np + tmpi
  end do
  temps%d(:) = 1
  n0 = 0.1d0
  print *, ''
  print *, 'Low-rank method (spherically symmetric speeds):'
  call MonteDSMCSphere(n0, nv, temps, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1)
!   call nv%deinit()
!   call temps%deinit()
!   call nv%init(2**10)
!   call temps%init(2**10)
!   nv%d(1) = 999
!   i = 1
!   tmpi = 1.0d0
!   np = nv%d(1)
!   do while (tmpi > 0)
!     i = i + 1
!     tmpi = floor(nv%d(1)*exp(-0.01d0*(i-1)))
!     nv%d(i) = floor(nv%d(1)*exp(-0.01d0*(i-1)))
!     np = np + tmpi
!   end do
!   temps%d(:) = 1
!   n0 = 0.1d0
!   print *, ''
!   print *, 'Acceptance-rejection method:'
!   call DSMCup(n0, nv, temps, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1) !ONLY FOR NO AGGREGATION
end
