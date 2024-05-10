Module ModMonte
  USE ModMtrx
  USE MerTwist, only : grnd, ggrnd, ggrnd2
  
  !Module for low-rank Monte-Carlo simulations
  !Contains MonteDSMC and MonteTemp
  
  !Function interfaces. See definitions below.
  abstract interface
    pure function UV_int(t, i, r) Result(res)
      Double precision, intent(in) :: t
      Integer(4), intent(in) :: i, r
      Double precision :: res
    end
  end interface
  abstract interface
    subroutine Kernel_int(n, u, v)
      USE ModMtrx
      Integer(4), intent(in) :: n
      Type(Mtrx), intent(out) :: u, v
    end
  end interface
  abstract interface
    pure function e_int(vx, vy, vz, dx, dy, dz, m1, m2) Result(e)
      Double precision, intent(in) :: vx, vy, vz, dx, dy, dz
      Integer(4), intent(in) :: m1, m2
      Double precision :: e
    end
  end interface
  abstract interface
    pure function bounce_int(vx, vy, vz, dx, dy, dz, m1, m2) Result(bounce)
      Double precision, intent(in) :: vx, vy, vz, dx, dy, dz
      Integer(4), intent(in) :: m1, m2
      Logical :: bounce
    end
  end interface
  abstract interface
    pure function rate_int(i, j, ti, tj) Result(res)
      Integer(4), intent(in) :: i, j
      Double precision, intent(in) :: ti, tj
      Double precision :: res
    end
  end interface
  abstract interface
    pure function prob_int(lambda, a, i, j, ti, tj) Result(res)
      Double precision, intent(in) :: lambda, a
      Integer(4), intent(in) :: i, j
      Double precision, intent(in) :: ti, tj
      Double precision :: res
    end
  end interface
  abstract interface
    pure function plus_int(lambda, a, i, j, ti, tj, tk, nk) Result(res)
      Double precision, intent(in) :: lambda, a
      Integer(4), intent(in) :: i, j
      Double precision, intent(in) :: ti, tj, tk
      Integer(4), intent(in) :: nk
      Double precision :: res
    end
  end interface
  abstract interface
    pure function minus_int(lambda, a, i, j, ti, tj, ni) Result(res)
      Double precision, intent(in) :: lambda, a
      Integer(4), intent(in) :: i, j
      Double precision, intent(in) :: ti, tj
      Integer(4), intent(in) :: ni
      Double precision :: res
    end
  end interface
  
  abstract interface
    pure function maxrate(imin, jmin, imax, jmax) Result(res)
      Integer(4), intent(in) :: imin, jmin, imax, jmax
      Double precision :: res
    end
  end interface
  abstract interface
    pure function rate(i, j) Result(res)
      Integer(4), intent(in) :: i, j
      Double precision :: res
    end
  end interface
  abstract interface
    subroutine rate_uv(i, ti, u, v)
      USE ModVec
      Integer(4), intent(in) :: i
      Double precision, intent(in) :: ti
      Type(Vector) :: u, v
    end
  end interface
  abstract interface
    pure function rate_vec(i) Result(res)
      Integer(4), intent(in) :: i
      Double precision :: res
    end
  end interface
  
  !Basic functions, can be assigned by user. See functions with 0 in name for description.
  !For DSMC
  public :: MonteDSMC, MonteDSMCsphere!, DSMCup
  procedure(UV_int), pointer :: Uspeed => Uspeed0, Vspeed => Vspeed0
  procedure(Kernel_int), pointer :: DSMCKernel => DSMCKernel0
  procedure(e_int), pointer :: ecoef => ecoef0
  procedure(bounce_int), pointer :: BounceIf => BounceIf0
  private :: FullApproximate, FullUpdate, OneUpdate
  !For temperature-dependent Monte Carlo (MonteTemp)
  public :: MonteTemp, MonteUpBtemp!, MonteWagnerRanked
  procedure(UV_int), pointer :: Utemp => Utemp0, Vtemp => Vtemp0
  procedure(Kernel_int), pointer :: TempKernel => TempKernel0
  procedure(minus_int), pointer :: AggloMinus => AggloMinus0, BounceMinus => BounceMinus0
  procedure(plus_int), pointer :: AggloPlus => AggloPlus0
  procedure(rate_int), pointer :: sucrate => sucrate0, CijTemp => CijTemp0
  procedure(prob_int), pointer :: BounceProb => BounceProb0
  private :: TempApproximate, TempUpdate
  !For Walker method
  public :: MonteWalker
  procedure(maxrate), pointer :: Cmax => Cmax0
  procedure(rate), pointer :: Cij => Cij0, sucrates => sucrates0
  private :: ApproximateNo, initWalker, getWalker
  !For Wagner method
  public :: MonteWagnerSimple
  procedure(rate_uv), pointer :: Cuv => Cuv0
  procedure(rate_vec), pointer :: u_vec => u_vec0, v_vec => v_vec0
  !Low-rank for classical Smoluchowski equations
  public :: MonteSimple, MonteSimplePart
  !Classical methods
  public :: MonteUpB, MonteInverse
  
  Contains
  
  !ASSIGNABLE SUBROUTINES. Examples with descriptions for interfaces above.
  
!Collision kernel
subroutine DSMCKernel0(n, u, v)
  Integer(4), intent(in) :: n !Kernel size
  Type(Mtrx), intent(out) :: u, v !Kernel low-rank factors
  
  Integer(4) i, r

  r = 3
  call u%init(n, r)
  call v%init(n, r)
  !u * v = 2 sigma_{ij} = 2 (i^{1/3}/2 + j^{1/3}/2)^2
  do i = 1, n
    v%d(i,1) = 0.5d0
    v%d(i,2) = i**(1.0d0/3.0d0)
    v%d(i,3) = i**(2.0d0/3.0d0)/2
    u%d(i,1) = i**(2.0d0/3.0d0)
    u%d(i,2) = i**(1.0d0/3.0d0)
    u%d(i,3) = 1.0d0
  end do
end

!Speed part of U factor
pure function Uspeed0(s, i, r) Result(res)
  Double precision, intent(in) :: s !Speed
  Integer(4), intent(in) :: i, r !Row and column
  Double precision :: res !Multiplier

  res = s
end

!Speed part of V factor
pure function Vspeed0(s, i, r) Result(res)
  Double precision, intent(in) :: s !Speed
  Integer(4), intent(in) :: i, r !Row and column
  Double precision :: res !Multiplier

  res = 1.0d0
end

!Restitution coefficient
pure function ecoef0(vx, vy, vz, dx, dy, dz, m1, m2) Result(e)
  Double precision, intent(in) :: vx, vy, vz, dx, dy, dz !Relative speed v and collision direction d
  Integer(4), intent(in) :: m1, m2 !Particle masses
  Double precision :: e !Restitution coefficient
  
  e = 0.99d0
  
  !e = 1.0d0
end

!Bouncing condition
pure function BounceIf0(vx, vy, vz, dx, dy, dz, m1, m2) Result(bounce)
  Double precision, intent(in) :: vx, vy, vz, dx, dy, dz !Relative speed v and collision direction d
  Integer(4), intent(in) :: m1, m2 !Particle masses
  Logical :: bounce !TRUE if bounce; FALSE if aggregation
  
  !Double precision meff, Ekin
  
  bounce = .true.
  
  !bounce = .false.
  
  !meff = 1.0d0/(1.0d0/m1+1.0d0/m2)
  !Ekin = meff/2*(vx**2+vy**2+vz**2)
  !if (Ekin*0.99d0**2 > (0.15d0*0.99d0**2)*(1.0d0/(1.0d0/m1**(1.0d0/3.0d0)+1.0d0/m2**(1.0d0/3.0d0)))**1.0d0) then
  !  bounce = .true.
  !else
  !  bounce = .false.
  !end if
end





!Collision kernel
subroutine TempKernel0(n, u, v)
  Integer(4), intent(in) :: n !Kernel size
  Type(Mtrx), intent(out) :: u, v !Kernel low-rank factors
  
  Integer(4) i, r

  r = 3
  call u%init(n, r)
  call v%init(n, r)
  !sigma_{ij} = (i^{1/3}/2 + j^{1/3}/2)^2/sqrt(i)
  !Also devided by sqrt(i) to account for sqrt(T_i/i)
  !sqrt(T_i) is accounted for later
  do i = 1, n
    v%d(i,1) = 1.0d0*sqrt(pi/2)
    v%d(i,2) = 2*i**(1.0d0/3.0d0)*sqrt(pi/2)
    v%d(i,3) = i**(2.0d0/3.0d0)*sqrt(pi/2)
    u%d(i,1) = i**(1.0d0/6.0d0)
    u%d(i,2) = i**(-1.0d0/6.0d0)
    u%d(i,3) = i**(-1.0d0/2.0d0)
  end do
end

!Temperature part of U factor
pure function Utemp0(t, i, r) Result(res)
  Double precision, intent(in) :: t !Temperature
  Integer(4), intent(in) :: i, r !Row and column
  Double precision :: res !Multiplier

  res = sqrt(t)
end

!Temperature part of V factor
pure function Vtemp0(t, i, r) Result(res)
  Double precision, intent(in) :: t !Temperature
  Integer(4), intent(in) :: i, r !Row and column
  Double precision :: res !Multiplier

  res = 1.0d0
end

!Success rate of Cij upper bound
pure function sucrate0(i, j, ti, tj) Result(res)
  Integer(4), intent(in) :: i, j !Masses
  Double precision, intent(in) :: ti, tj !Temperatures
  Double precision :: res !Success probability
  
  res = sqrt(ti/i + tj/j)/(sqrt(ti/i) + sqrt(tj/j))
end

!Temperature-dependent kernel values as-is
pure function CijTemp0(i,j,ti,tj) Result(res)
  Integer(4), intent(in) :: i, j !Masses
  Double precision, intent(in) :: ti, tj !Temperatures
  Double precision :: res !Collision rate
  
  res = (i**(1.0d0/3.0d0)+j**(1.0d0/3.0d0))**2*sqrt(ti/i+tj/j)
end

!Probability of bouncing (restitution)
pure function BounceProb0(lambda, a, i, j, ti, tj) Result(res)
  Double precision, intent(in) :: lambda, a !Parameters
  Integer(4), intent(in) :: i, j !Masses
  Double precision, intent(in) :: ti, tj !Temperatures
  Double precision :: res !Resulting probability
  
  Double precision e, lambda1, lambda2, i3, j3, w, q
  
  if ((i .ne. 0) .and. (j .ne. 0) .and. (ti + tj .ne. 0.0d0)) then
    e = 0.99d0
    lambda1 = lambda
    lambda2 = lambda
    i3 = i**(1.0d0/3.0d0)
    j3 = j**(1.0d0/3.0d0)
    w = a * (i3*j3)**lambda1 / (i3 + j3)**lambda2
    q = w/e**2*(i+j)/(ti*j+tj*i)
  
    res = exp(-q)*(1+q)
  else
    res = 1.0d0
  end if
end

!Temperature change of k = i + j particle
pure function AggloPlus0(lambda, a, i, j, ti, tj, tk, nk) Result(res)
  Double precision, intent(in) :: lambda, a !Parameters
  Integer(4), intent(in) :: i, j !Masses
  Double precision, intent(in) :: ti, tj, tk !Temperatures
  Integer(4), intent(in) :: nk !Number of particles of size k
  Double precision :: res !New temperature
  
  Double precision e, lambda1, lambda2, i3, j3, w, q, f, g
  e = 0.99d0
  lambda1 = lambda
  lambda2 = lambda
  i3 = i**(1.0d0/3.0d0)
  j3 = j**(1.0d0/3.0d0)
  w = a * (i3*j3)**lambda1 / (i3 + j3)**lambda2
  q = w/e**2*(i+j)/(ti*j+tj*i)
  f = exp(-q)*(1+q)
  g = exp(-q)*(1+q+q**2/2)
  
  res = (nk*tk + (ti*tj*(1.0d0/i+1.0d0/j) + 4.0d0/3.0d0*(ti-tj)**2*(1-g)/(1-f)/(i+j))/(ti/i + tj/j))/(nk + 1)
end

!Temperature change of size i particles, when they aggregate
pure function AggloMinus0(lambda, a, i, j, ti, tj, ni) Result(res)
  Double precision, intent(in) :: lambda, a !Parameters
  Integer(4), intent(in) :: i, j !Masses
  Double precision, intent(in) :: ti, tj !Temperatures
  Integer(4), intent(in) :: ni !Number of size i particles
  Double precision :: res !New temperature
  
  Double precision e, lambda1, lambda2, i3, j3, w, q, f, g
  e = 0.99d0
  lambda1 = lambda
  lambda2 = lambda
  i3 = i**(1.0d0/3.0d0)
  j3 = j**(1.0d0/3.0d0)
  w = a * (i3*j3)**lambda1 / (i3 + j3)**lambda2
  q = w/e**2*(i+j)/(ti*j+tj*i)
  f = exp(-q)*(1+q)
  g = exp(-q)*(1+q+q**2/2)
  if (ni > 1) then
    res = ti*(ni - (tj/j + 4.0d0/3.0d0*ti/i*(1-g)/(1-f))/(ti/i+tj/j))/(ni - 1)
    res = max(res, 0.0d0)
  else
    res = 0.0d0
  end if
end

!Temperature change of size i particles, when they bounce
pure function BounceMinus0(lambda, a, i, j, ti, tj, ni) Result(res)
  Double precision, intent(in) :: lambda, a !Parameters
  Integer(4), intent(in) :: i, j !Masses
  Double precision, intent(in) :: ti, tj !Temperatures
  Integer(4), intent(in) :: ni !Number of size i particles
  Double precision :: res !New temperature
  
  Double precision e, lambda1, lambda2, i3, j3, w, q, f, g
  e = 0.99d0
  lambda1 = lambda
  lambda2 = lambda
  i3 = i**(1.0d0/3.0d0)
  j3 = j**(1.0d0/3.0d0)
  w = a * (i3*j3)**lambda1 / (i3 + j3)**lambda2
  q = w/e**2*(i+j)/(ti*j+tj*i)
  f = exp(-q)*(1+q)
  g = exp(-q)*(1+q+q**2/2)

  res = 4.0d0/3.0d0*j*(1+e)*(ti*(i+j) - 0.5d0*(1+e)*(ti*j + tj*i))/dble(i+j)**2/ni
  res = res*g/f
  res = ti - res
  res = max(res, 0.0d0)
end





!Upper bound of Cij
pure function Cmax0(imin, jmin, imax, jmax) Result(res)
  Integer(4), intent(in) :: imin, jmin, imax, jmax !Rectangle, where maximum is searched
  Double precision :: res !Upper bound on collision rates
  
  !Linear kernel
  res = dble(imax) + jmax
  !Brownian kernel
  !res = max(dble(imax)/jmin,dble(jmax)/imin)
  !res = 2+res**(1.0d0/3.0d0)+res**(-1.0d0/3.0d0)
  !Ballistic kernel
  !res = (imin**(1.0d0/3.0d0) + jmax**(1.0d0/3.0d0))**2*sqrt(1.0d0/imin+1.0d0/jmax)
  !res = max(res, (imax**(1.0d0/3.0d0) + jmin**(1.0d0/3.0d0))**2*sqrt(1.0d0/imax+1.0d0/jmin))
end

!Classical kernel Cij
pure function Cij0(i, j) Result(res)
  Integer(4), intent(in) :: i, j !Masses
  Double precision :: res !Collision rate
  
  res = dble(i) + j
end

!Success rate of Cij upper bound
pure function sucrates0(i, j) Result(res)
  Integer(4), intent(in) :: i, j !Masses
  Double precision :: res !Success probability
  
  res = 1.0d0
end






subroutine Cuv0(i,ti,u,v)
  USE ModVec
  Integer(4), intent(in) :: i
  Double precision, intent(in) :: ti
  Type(Vector) :: u, v

  u%d(1) = i**(2.0d0/3.0d0)
  u%d(2) = 2.0d0*i**(1.0d0/3.0d0)
  u%d(3) = 1.0d0
  v%d(1) = sqrt(ti/i)
  v%d(2) = i**(1.0d0/3.0d0)*sqrt(ti/i)
  v%d(3) = i**(2.0d0/3.0d0)*sqrt(ti/i)
end

pure function u_vec0(i) Result(res)
  Integer(4), intent(in) :: i
  Double precision :: res
  
  res = i
  !res = 2*i**(1.0d0/3.0d0)
  !res = i**(2.0d0/3.0d0)
  !res = i**(0.98d0)
end

pure function v_vec0(i) Result(res)
  Integer(4), intent(in) :: i
  Double precision :: res
  
  res = 1.0d0
  !res = i**(-1.0d0/3.0d0)
  !res = 2*sqrt(2.0d0/i)
  !res = i**(-0.98d0)
end
  
  !GENERAL SUBROUTINES

!General DSMC calculation with ballistic trajectories.
subroutine MonteDSMC(n0, partofsize, temps, maxtime, gbig, gsmall, Tmolgas, verbose)
  Double precision, intent(inout) :: n0 !Total number density; IN and OUT
  Type(IntVec) :: partofsize !Numbers of particles; IN and OUT
  Type(Vector) :: temps !Temperatures of each size; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Double precision, intent(in) :: gbig, gsmall, tmolgas !Thermostat parameters
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every Nh steps
  
  Integer(4) np, np0 !Current and initial total number of particles
  Integer(4), parameter :: partmin = 4 !Minimum array size for particles of each mass
  Double precision e !Restitution coefficient
  Double precision tavg !Current average temperature
  Type(Vector) vx, vy, vz !Particle speeds
  Type(Vector) vs !Particle squared speed
  Type(Vector) vmaxvec !Maximum speed for each size
  Type(Mtrx) u, v, ui, vi !Low-rank factors
  Integer(4), allocatable :: sizestart(:) !Start index of size-i particles in the initial array
  
  Double precision curtime !Current system time
  Double precision ex, ey, ez !Collision direction components
  Double precision dx, dy, dz, vcs, vcs1, vcs2 !Temporaries for speed values
  Double precision cursec !Current time
  Double precision dsecnd !LAPACK procedure for time measurement
  Integer(4) maxmass !Maximum cluster mass in the system
  Integer(4) treestart !Start index of segment tree data
  Double precision ran, ran1, ran2, randi !Random numbers
  Double precision totalrate !Total collision rate
  Integer(4) r, curr !Rank, iteration over ranks
  Type(Vector) usum, vsum, totalvec !Vectors, containing sums of collision rates
  Type(Mtrx) utree, vtree !Segment trees
  Integer(8) Nh !Update output every Nh steps
  Integer(8) Nht !Apply thermostat every Nht steps
  Integer(8) step !Current step
  Double precision lasttime !System time at last thermostat application
  Double precision tempsum !Total kinetic energy (sum of weighted temperatures)
  Type(IntVec), allocatable :: partofnum(:) !Arrays of particles of each mass
  Type(IntVec) vmass !Array of particle masses
  Type(Vector) vmassvec !Real array of particle masses
  Type(Mtrx) sm !Helps to select collision direction
  Double precision timeold, timecur !Saves
  Integer(4) i, j, i0, j0, i1, j1, tmpi, tmpi2, i2, j2, k2 !Indices
  
  !Temporaries
  Double precision tau, tau2, tauexp, tmp
  Integer(4), allocatable :: tmpint(:)
  Type(IntVec), allocatable :: tmpintvec(:)
  Type(Vector) tmpvec, uioi, uioj, uiok, vioi, vioj, viok
  
  maxmass = max(16, ibset(0, 2*((bit_size(partofsize%n-1) - leadz(partofsize%n-1) + 1)/2)))
  Allocate(tmpint(maxmass))
  tmpint(1:partofsize%n) = partofsize%d(1:partofsize%n)
  tmpint(partofsize%n+1:) = 0
  partofsize = tmpint
  Deallocate(tmpint)
  call tmpvec%copy(temps)
  call temps%deinit()
  call temps%init(maxmass)
  temps%d(1:tmpvec%n) = tmpvec%d(1:tmpvec%n)
  call tmpvec%deinit()
  where(partofsize%d==0) temps%d = 0
  
  np = sum(partofsize%d)
  np0 = np
  
  !Update time
  Nh = np!/20!5
  !Thermostat time
  Nht = np!/100
  
  Allocate(sizestart(maxmass+1))
  Allocate(partofnum(maxmass))
  call vmass%init(np+1)
  call vmassvec%init(np)
  sizestart(1) = 1
  i = 1
  j1 = 1
  do while (i <= maxmass)
    tmpi = partofsize%d(i)
    call partofnum(i)%init(max(tmpi,partmin))
    do j = j1, j1+tmpi-1
      vmass%d(j) = i
      partofnum(i)%d(j-j1+1) = j
    end do
    j1 = j1 + tmpi
    i = i + 1
    sizestart(i) = sizestart(i-1) + tmpi
  end do
  
  treestart = 0
  !Calculate geometric series
  i1 = 1
  do while (i1*4 < maxmass)
    treestart = treestart + i1
    i1 = i1*4
  end do

  call vx%init(np)
  call vy%init(np)
  call vz%init(np)
  call vs%init(np+1)
  call sm%init(3,2)
  
  do i = 1, np
    vx%d(i) = ggrnd()
    vy%d(i) = ggrnd()
    vz%d(i) = ggrnd()
  end do
  
  do j = 1, maxmass
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      ex = ex + vx%d(i)
      ey = ey + vy%d(i)
      ez = ez + vz%d(i)
    end do
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vx%d(i) = vx%d(i) - ex/partofsize%d(j)
      vy%d(i) = vy%d(i) - ey/partofsize%d(j)
      vz%d(i) = vz%d(i) - ez/partofsize%d(j)
    end do
    tmp = 0.0d0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      tmp = tmp + vx%d(i)**2 + vy%d(i)**2 + vz%d(i)**2
    end do
    tmp = tmp*j
    if (partofsize%d(j) > 1) then
      tmp = 3*temps%d(j)/tmp*partofsize%d(j)
    elseif (partofsize%d(j) == 1) then
      call randir(vx%d(sizestart(j)), vy%d(sizestart(j)), vz%d(sizestart(j)))
      tmp = 3*temps%d(j)/j
    end if
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vx%d(i) = vx%d(i)*sqrt(tmp)
      vy%d(i) = vy%d(i)*sqrt(tmp)
      vz%d(i) = vz%d(i)*sqrt(tmp)
    end do
  end do
  
  do i = 1, np
    vx%d(i) = vx%d(i)/sqrt(temps%d(vmass%d(i)))
    vy%d(i) = vy%d(i)/sqrt(temps%d(vmass%d(i)))
    vz%d(i) = vz%d(i)/sqrt(temps%d(vmass%d(i)))
  end do
  do i = 1, np
    vs%d(i) = vx%d(i)**2+vy%d(i)**2+vz%d(i)**2
  end do
  tempsum = 0
  do j = 1, maxmass
    temps%d(j) = 0.0d0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      temps%d(j) = temps%d(j) + vs%d(i)
    end do
    tempsum = tempsum + temps%d(j)*j
    temps%d(j) = (temps%d(j)*j)/max(1,partofsize%d(j))/3
  end do
  if (verbose > 0) then
    print *, 'Monomer Temperature', temps%d(1)
    print *, 'TOTAL ENERGY START', tempsum
  end if
  
  call vmaxvec%init(maxmass)
  do j = 1, maxmass
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vmaxvec%d(j) = max(vmaxvec%d(j), vs%d(i))
    end do
    vmaxvec%d(j) = sqrt(vmaxvec%d(j))
  end do
  
  call FullApproximate(u, v, ui, vi, partofsize, vmaxvec)
  r = u%m
  call usum%init(r)
  call vsum%init(r)
  call uioi%init(r)
  call uioj%init(r)
  call uiok%init(r)
  call vioi%init(r)
  call vioj%init(r)
  call viok%init(r)
  call totalvec%init(r)
  utree = treebuild(ui, treestart)
  vtree = treebuild(vi, treestart)
  usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
  vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
  totalvec%d(:) = usum%d(:)*vsum%d(:)
  totalrate = sum(totalvec%d)
  
  curtime = 0.0d0
  lasttime = 0.0d0
  
  timeold = curtime
  
  step = 0
!   open(2, file='TimeTemps3')
!   open(3, file='TimeConc3')
  cursec = dsecnd()
  do while (curtime < maxtime)
    step = step + 1
    i1 = 0
    do while (1 > 0)
      !2 comes from integration over e vector in advance
      curtime = curtime + 1.0d0/n0*np0/totalrate/pi*2
      
      randi = grnd()*totalrate
      do curr = 1, r-1
        randi = randi - totalvec%d(curr)
        if (randi <= 0) then
          exit
        end if
      end do
      randi = grnd()*usum%d(curr)
      i2 = treefindr(utree, ui, randi, curr, treestart)
      randi = grnd()*vsum%d(curr)
      j2 = treefindr(vtree, vi, randi, curr, treestart)
      
      tmpi = floor(grnd()*partofsize%d(i2))+1
      tmpi2 = floor(grnd()*partofsize%d(j2))+1
      if (tmpi2 > tmpi) then
        i0 = tmpi
        j0 = tmpi2
      else
        i0 = tmpi2
        j0 = tmpi
        tmpi2 = i2
        i2 = j2
        j2 = tmpi2
      end if
      
      if ((j0 <= partofsize%d(j2)) .and. (i0 <= partofsize%d(i2))) then
        i1 = partofnum(i2)%d(i0)
        j1 = partofnum(j2)%d(j0)
        dx = vx%d(i1)-vx%d(j1)
        dy = vy%d(i1)-vy%d(j1)
        dz = vz%d(i1)-vz%d(j1)
        vcs = dx**2 + dy**2 + dz**2
        if (vcs > (grnd()*(vmaxvec%d(i2)+vmaxvec%d(j2)))**2) then
          exit
        end if
      end if
    end do
    
    vcs1 = sqrt(vcs)
    if (abs(dy) > abs(dx)) then
      vcs2 = sqrt(dz**2 + dy**2)
      sm%d(1,1) = 0
      sm%d(2,1) = -dz/vcs2
      sm%d(3,1) = dy/vcs2
      sm%d(1,2) = vcs2/vcs1
      sm%d(2,2) = -dx*dy/vcs2/vcs1
      sm%d(3,2) = -dx*dz/vcs2/vcs1
    else
      vcs2 = sqrt(dz**2 + dx**2)
      sm%d(1,1) = -dz/vcs2
      sm%d(2,1) = 0
      sm%d(3,1) = dx/vcs2
      sm%d(1,2) = -dx*dy/vcs2/vcs1
      sm%d(2,2) = vcs2/vcs1
      sm%d(3,2) = -dz*dy/vcs2/vcs1
    end if
    ran = sqrt(grnd())
    ran1 = sqrt(1 - ran**2)
    ran2 = grnd()*2*pi
    tau = ran1*cos(ran2)
    tau2 = ran1*sin(ran2)
    ex = ran*dx/vcs1 + tau*sm%d(1,1) + tau2*sm%d(1,2)
    ey = ran*dy/vcs1 + tau*sm%d(2,1) + tau2*sm%d(2,2)
    ez = ran*dz/vcs1 + tau*sm%d(3,1) + tau2*sm%d(3,2)
    
    if (BounceIf(dx, dy, dz, ex, ey, ez, i2, j2)) then
      e = ecoef(dx, dy, dz, ex, ey, ez, i2, j2)
      vcs = (1+e)*vcs1*ran
      vcs1 = vcs*j2/(i2 + j2)
      vcs2 = vcs*i2/(i2 + j2)
      vx%d(i1) = vx%d(i1) - vcs1*ex
      vy%d(i1) = vy%d(i1) - vcs1*ey
      vz%d(i1) = vz%d(i1) - vcs1*ez
      vx%d(j1) = vx%d(j1) + vcs2*ex
      vy%d(j1) = vy%d(j1) + vcs2*ey
      vz%d(j1) = vz%d(j1) + vcs2*ez
      vcs = vx%d(i1)**2+vy%d(i1)**2+vz%d(i1)**2
      vs%d(i1) = vcs
      if (vcs > vmaxvec%d(i2)**2) then
        vmaxvec%d(i2) = sqrt(vcs)
        
        uioi%d(:) = ui%d(i2,:)
        vioi%d(:) = vi%d(i2,:)
        
        call OneUpdate(i2, u, v, ui, vi, partofsize, vmaxvec)
        
        uioi%d(:) = ui%d(i2,:) - uioi%d(:)
        vioi%d(:) = vi%d(i2,:) - vioi%d(:)
        
        call treeupdate2(utree, ui, i2, treestart, uioi)
        call treeupdate2(vtree, vi, i2, treestart, vioi)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
      end if
      vcs = vx%d(j1)**2+vy%d(j1)**2+vz%d(j1)**2
      vs%d(j1) = vcs
      if (vcs > vmaxvec%d(j2)**2) then
        vmaxvec%d(j2) = sqrt(vcs)
        
        uioi%d(:) = ui%d(j2,:)
        vioi%d(:) = vi%d(j2,:)
        
        call OneUpdate(j2, u, v, ui, vi, partofsize, vmaxvec)
        
        uioi%d(:) = ui%d(j2,:) - uioi%d(:)
        vioi%d(:) = vi%d(j2,:) - vioi%d(:)
        
        call treeupdate2(utree, ui, j2, treestart, uioi)
        call treeupdate2(vtree, vi, j2, treestart, vioi)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
      end if
    else
      partofnum(j2)%d(j0) = partofnum(j2)%d(partofsize%d(j2))
      partofsize%d(j2) = partofsize%d(j2) - 1
      partofnum(i2)%d(i0) = partofnum(i2)%d(partofsize%d(i2))
      partofsize%d(i2) = partofsize%d(i2) - 1
      k2 = i2 + j2
      if (k2 > maxmass) then
        Allocate(tmpintvec(maxmass))
        tmpintvec(:) = partofnum(:)
        Deallocate(partofnum)
        Allocate(partofnum(4*maxmass))
        partofnum(1:maxmass) = tmpintvec(1:maxmass)
        do j = 1+maxmass,4*maxmass
          call partofnum(j)%init(partmin)
        end do
        Deallocate(tmpintvec)
        Allocate(tmpint(maxmass))
        tmpint(:) = partofsize%d(:)
        call partofsize%deinit()
        call partofsize%init(4*maxmass)
        partofsize%d(1:maxmass) = tmpint(1:maxmass)
        Deallocate(tmpint)
        call tmpvec%init(maxmass)
        tmpvec%d(:) = vmaxvec%d(:)
        call vmaxvec%deinit()
        call vmaxvec%init(4*maxmass)
        vmaxvec%d(1:maxmass) = tmpvec%d(1:maxmass)
        call tmpvec%deinit()
        maxmass = maxmass*4
        
        call u%deinit()
        call v%deinit()
        call ui%deinit()
        call vi%deinit()
        call FullApproximate(u, v, ui, vi, partofsize, vmaxvec)
        treestart = treestart + maxmass/4
        utree = treebuild(ui, treestart)
        vtree = treebuild(vi, treestart)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
        
      end if
      if (partofsize%d(k2) >= partofnum(k2)%n) then
        Allocate(tmpint(partofsize%d(k2)))
        tmpint(:) = partofnum(k2)%d(:)
        call partofnum(k2)%deinit()
        call partofnum(k2)%init(2*partofsize%d(k2))
        partofnum(k2)%d(1:partofsize%d(k2)) = tmpint(1:partofsize%d(k2))
        Deallocate(tmpint)
      end if
      partofsize%d(k2) = partofsize%d(k2) + 1
      partofnum(k2)%d(partofsize%d(k2)) = i1
      
      vx%d(i1) = (i2*vx%d(i1)+j2*vx%d(j1))/k2
      vy%d(i1) = (i2*vy%d(i1)+j2*vy%d(j1))/k2
      vz%d(i1) = (i2*vz%d(i1)+j2*vz%d(j1))/k2
      vmass%d(i1) = k2
      vmass%d(j1) = 0
      vcs = vx%d(i1)**2+vy%d(i1)**2+vz%d(i1)**2
      vs%d(i1) = vcs
      if (vcs > vmaxvec%d(k2)**2) then
        vmaxvec%d(k2) = sqrt(vcs)
      end if
      
      uioi%d(:) = ui%d(i2,:)
      vioi%d(:) = vi%d(i2,:)
      
      call OneUpdate(i2, u, v, ui, vi, partofsize, vmaxvec)
      
      uioi%d(:) = ui%d(i2,:) - uioi%d(:)
      vioi%d(:) = vi%d(i2,:) - vioi%d(:)
      
      call treeupdate2(utree, ui, i2, treestart, uioi)
      call treeupdate2(vtree, vi, i2, treestart, vioi)
      
      uioj%d(:) = ui%d(j2,:)
      vioj%d(:) = vi%d(j2,:)
      
      call OneUpdate(j2, u, v, ui, vi, partofsize, vmaxvec)
      
      uioj%d(:) = ui%d(j2,:) - uioj%d(:)
      vioj%d(:) = vi%d(j2,:) - vioj%d(:)
      
      call treeupdate2(utree, ui, j2, treestart, uioj)
      call treeupdate2(vtree, vi, j2, treestart, vioj)
      
      uiok%d(:) = ui%d(k2,:)
      viok%d(:) = vi%d(k2,:)
      
      call OneUpdate(k2, u, v, ui, vi, partofsize, vmaxvec)
      
      uiok%d(:) = ui%d(k2,:) - uiok%d(:)
      viok%d(:) = vi%d(k2,:) - viok%d(:)
      
      call treeupdate2(utree, ui, k2, treestart, uiok)
      call treeupdate2(vtree, vi, k2, treestart, viok)
      
      usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
      vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
      totalvec%d(:) = usum%d(:)*vsum%d(:)
      totalrate = sum(totalvec%d)
      
      np = np - 1
      
      !Number of particles doubling
      if (np <= np0/2) then
        np = np*2
        i = 1
        do j = 1, maxmass
          j1 = partofnum(j)%n
          if (2*partofsize%d(j) > j1) then
            Allocate(tmpint(partofsize%d(j)))
            tmpint(1:partofsize%d(j)) = partofnum(j)%d(1:partofsize%d(j))
            call partofnum(j)%deinit()
            call partofnum(j)%init(2*j1)
            partofnum(j)%d(1:partofsize%d(j)) = tmpint(1:partofsize%d(j))
            Deallocate(tmpint)
          end if
          do j1 = 1 + partofsize%d(j), 2*partofsize%d(j)
            do while(vmass%d(i) .ne. 0)
              i = i + 1
            end do
            vmass%d(i) = j
            j2 = partofnum(j)%d(j1-partofsize%d(j))
            call randir(ex, ey, ez)
            vs%d(i) = vs%d(j2)
            vx%d(i) = ex*sqrt(vs%d(j2))
            vy%d(i) = ey*sqrt(vs%d(j2))
            vz%d(i) = ez*sqrt(vs%d(j2))
            partofnum(j)%d(j1) = i
          end do
        end do
        partofsize%d(:) = partofsize%d(:)*2
        ui = ui*2.0d0
        vi = vi*2.0d0
        utree = treebuild(ui, treestart)
        vtree = treebuild(vi, treestart)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
        n0 = n0/2
      end if
    end if
    
    !Thermostat (molecular gas)
    if ((mod(step,Nht) == 0) .and. (gbig > 0)) then
      tau = gbig*(curtime-lasttime)
      do j2 = 1, maxmass
        if (Tmolgas > 0) then
          tau2 = sqrt(3*Tmolgas/j2*(1 - exp(-tau*(j2**(gsmall+1))/Tmolgas)))
          tauexp = exp(-tau*(j2**(gsmall+1))/Tmolgas/2)
        else
          tau2 = sqrt(3*tau*j2**(gsmall))
          tauexp = 1.0d0
        end if
        do j1 = 1, partofsize%d(j2)
          j = partofnum(j2)%d(j1)
          call randiradv(ex, ey, ez)
          vx%d(j) = vx%d(j)*tauexp+ex*tau2
          vy%d(j) = vy%d(j)*tauexp+ey*tau2
          vz%d(j) = vz%d(j)*tauexp+ez*tau2
        end do
      end do
      lasttime = curtime
    
      if (mod(step,Nh) .ne. 0) then
        vmaxvec%d = 0
        do j = 1, np0
          j1 = vmass%d(j)
          if (j1 > 0) then
            vmaxvec%d(j1) = max(vmaxvec%d(j1), vs%d(j))
          end if
        end do
        vmaxvec%d(1:maxmass) = sqrt(vmaxvec%d(1:maxmass))
      
        call FullUpdate(u, v, ui, vi, partofsize, vmaxvec)
        utree = treebuild(ui, treestart)
        vtree = treebuild(vi, treestart)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
      end if
    end if
    
    if ((mod(step,Nh) == 0) .or. (curtime >= maxtime)) then
      vmaxvec%d = 0
      tempsum = 0
      call temps%deinit()
      call temps%init(maxmass)
      do j = 1, np0
        j1 = vmass%d(j)
        if (j1 > 0) then
          vmaxvec%d(j1) = max(vmaxvec%d(j1), vs%d(j))
          temps%d(j1) = temps%d(j1) + vs%d(j)*j1
          tempsum = tempsum + vs%d(j)*j1
        end if
      end do
      vmaxvec%d(1:maxmass) = sqrt(vmaxvec%d(1:maxmass))
      temps%d(1:maxmass) = temps%d(1:maxmass)/max(1,partofsize%d(1:maxmass))/3
      if (verbose > 1) then
        print *, 'TOTAL ENERGY', tempsum
      end if
      tavg = tempsum/np/3
      
      call FullUpdate(u, v, ui, vi, partofsize, vmaxvec)
      utree = treebuild(ui, treestart)
      vtree = treebuild(vi, treestart)
      usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
      vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
      totalvec%d(:) = usum%d(:)*vsum%d(:)
      totalrate = sum(totalvec%d)

      if (verbose > 1) then
!         write(2, *) curtime, temps%d(1), temps%d(3), temps%d(5), temps%d(8), tavg
!         write(3, *) curtime, dble(partofsize%d(1))/np0, (n0*np)/np0
      
        print *, step, curtime, np, maxmass
        print *, 'Temperature', temps%d(1), temps%d(2), temps%d(5), temps%d(8), tavg
        print *, 'Time passed', dsecnd()-cursec
        
        timecur = curtime
        timeold = timecur
        
      end if
    end if
  end do
!   close(2)
!   close(3)
  if (verbose > 0) then
    if (verbose == 1) then
      print *, 'Time passed', dsecnd()-cursec
    end if
    print *, 'Collisions', step
    n0 = (n0*np)/np0
    print *, 'Number density', n0
    print *, 'DONE'
  end if
  
  n0 = (n0*np)/np0
  
  Deallocate(sizestart)
  do j = 1, size(partofnum)
    call partofnum(j)%deinit()
  end do
  Deallocate(partofnum)
end

!DSMC with spherically symmetric speed destribution.
subroutine MonteDSMCsphere(n0, partofsize, temps, maxtime, gbig, gsmall, Tmolgas, verbose)
  Double precision, intent(inout) :: n0 !Total number density; IN and OUT
  Type(IntVec) :: partofsize !Numbers of particles; IN and OUT
  Type(Vector) :: temps !Temperatures of each size; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Double precision, intent(in) :: gbig, gsmall, tmolgas !Thermostat parameters
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every Nh steps
  
  Integer(4) np, np0 !Current and initial total number of particles
  Integer(4), parameter :: partmin = 4 !Minimum array size for particles of each mass
  Double precision e !Restitution coefficient
  Double precision tavg !Current average temperature
  Type(Vector) vx, vy, vz !Particle speeds
  Type(Vector) vs !Particle squared speed
  Type(Vector) vmaxvec !Maximum speed for each size
  Type(Mtrx) u, v, ui, vi !Low-rank factors
  Integer(4), allocatable :: sizestart(:) !Start index of size-i particles in the initial array
  
  Double precision curtime !Current system time
  Double precision ex, ey, ez !Collision direction components
  Double precision dx, vcs, vcs1, vcs2, vr, va, vb, v1, v2 !Temporaries for speed values
  Double precision cursec !Current time
  Double precision dsecnd !LAPACK procedure for time measurement
  Integer(4) maxmass !Maximum cluster mass in the system
  Integer(4) treestart !Start index of segment tree data
  Double precision ran, ran1, ran2, randi !Random numbers
  Double precision totalrate !Total collision rate
  Integer(4) r, curr !Rank, iteration over ranks
  Type(Vector) usum, vsum, totalvec !Vectors, containing sums of collision rates
  Type(Mtrx) utree, vtree !Segment trees
  Integer(8) Nh !Update output every Nh steps
  Integer(8) Nht !Apply thermostat every Nht steps
  Integer(8) step !Current step
  Double precision lasttime !System time at last thermostat application
  Double precision tempsum !Total kinetic energy (sum of weighted temperatures)
  Type(IntVec), allocatable :: partofnum(:) !Arrays of particles of each mass
  Type(IntVec) vmass !Array of particle masses
  Type(Vector) vmassvec !Real array of particle masses
  Double precision timeold, timecur !Saves
  Integer(4) i, j, i0, j0, i1, j1, tmpi, tmpi2, i2, j2, k2 !Indices
  Type(IntVec) ones !Vector of ones, to avoid division by zero.
  
  !Temporaries
  Double precision tau, tau2, tauexp, tmp, tmp2
  Integer(4), allocatable :: tmpint(:)
  Type(IntVec), allocatable :: tmpintvec(:)
  Type(Vector) tmpvec, uioi, uioj, uiok, vioi, vioj, viok
  
  maxmass = max(16, ibset(0, 2*((bit_size(partofsize%n-1) - leadz(partofsize%n-1) + 1)/2)))
  Allocate(tmpint(maxmass))
  tmpint(1:partofsize%n) = partofsize%d(1:partofsize%n)
  tmpint(partofsize%n+1:) = 0
  partofsize = tmpint
  Deallocate(tmpint)
  call tmpvec%copy(temps)
  call temps%deinit()
  call temps%init(maxmass)
  temps%d(1:tmpvec%n) = tmpvec%d(1:tmpvec%n)
  call tmpvec%deinit()
  where(partofsize%d == 0) temps%d = 0
  
  np = sum(partofsize%d)
  np0 = np
  tavg = temps%d(1)
  
  !Update time
  Nh = np!/20!5
  !Thermostat time
  Nht = np!/100
  
  Allocate(sizestart(maxmass+1))
  Allocate(partofnum(maxmass))
  call vmass%init(np+1)
  call vmassvec%init(np)
  sizestart(1) = 1
  i = 1
  j1 = 1
  do while (i <= maxmass)
    tmpi = partofsize%d(i)
    call partofnum(i)%init(max(tmpi,partmin))
    do j = j1, j1+tmpi-1
      vmass%d(j) = i
      partofnum(i)%d(j-j1+1) = j
    end do
    j1 = j1 + tmpi
    i = i + 1
    sizestart(i) = sizestart(i-1) + tmpi
  end do
  
  treestart = 0
  !Calculate geometric series
  i1 = 1
  do while (i1*4 < maxmass)
    treestart = treestart + i1
    i1 = i1*4
  end do

  call vx%init(np)
  call vy%init(np)
  call vz%init(np)
  call vs%init(np+1)
  do i = 1, np
    vs%d(i) = sqrt(ggrnd()**2 - log(grnd()))
  end do
  
  do j = 1, maxmass
    tmp = 0.0d0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      tmp = tmp + vs%d(i)**2
    end do
    tmp = tmp*j
    if (partofsize%d(j) >= 1) then
      tmp = 3*temps%d(j)/tmp*partofsize%d(j)
    end if
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vs%d(i) = vs%d(i)*sqrt(tmp)
    end do
  end do

  tempsum = 0
  do j = 1, maxmass
    temps%d(j) = 0.0d0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      temps%d(j) = temps%d(j) + vs%d(i)**2
    end do
    tempsum = tempsum + temps%d(j)*j
    temps%d(j) = (temps%d(j)*j)/max(1,partofsize%d(j))/3
  end do
  if (verbose > 0) then
    print *, 'Monomer Temperature', temps%d(1)
    print *, 'TOTAL ENERGY START', tempsum
  end if
  
  call vmaxvec%init(maxmass)
  do j = 1, maxmass
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vmaxvec%d(j) = max(vmaxvec%d(j), vs%d(i))
    end do
  end do
  
  call FullApproximate(u, v, ui, vi, partofsize, vmaxvec)
  r = u%m
  call usum%init(r)
  call vsum%init(r)
  call uioi%init(r)
  call uioj%init(r)
  call uiok%init(r)
  call vioi%init(r)
  call vioj%init(r)
  call viok%init(r)
  call totalvec%init(r)
  utree = treebuild(ui, treestart)
  vtree = treebuild(vi, treestart)
  usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
  vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
  totalvec%d(:) = usum%d(:)*vsum%d(:)
  totalrate = sum(totalvec%d)
  
  curtime = 0.0d0
  lasttime = 0.0d0
  
  timeold = curtime
  
  step = 0
!   open(2, file='TimeTempsSph')
!   open(3, file='TimeConcSph')
  cursec = dsecnd()
  do while (curtime < maxtime)
    step = step + 1
    i1 = 0
    do while (1 > 0)
      !2 comes from integration over e vector in advance
      curtime = curtime + 1.0d0/n0*np0/totalrate/pi*2
      
      randi = grnd()*totalrate
      do curr = 1, r-1
        randi = randi - totalvec%d(curr)
        if (randi <= 0) then
          exit
        end if
      end do
      randi = grnd()*usum%d(curr)
      i2 = treefindr(utree, ui, randi, curr, treestart)
      randi = grnd()*vsum%d(curr)
      j2 = treefindr(vtree, vi, randi, curr, treestart)
      
      i0 = floor(grnd()*partofsize%d(i2))+1
      j0 = floor(grnd()*partofsize%d(j2))+1
      
      if ((j0 <= partofsize%d(j2)) .and. (i0 <= partofsize%d(i2))) then
        i1 = partofnum(i2)%d(i0)
        j1 = partofnum(j2)%d(j0)
        if (i1 .ne. j1) then
          v1 = vs%d(i1)
          v2 = vs%d(j1)
          vr = v1 + v2 + abs(v1 - v2)
          if (vr**2 + 2*v1**2 + 2*v2**2 > grnd()*(vmaxvec%d(i2)+vmaxvec%d(j2))*3*vr) then
            exit
          end if
        end if
      end if
    end do
    
    va = (v1**2 + v2**2 - ((v1+v2)**3 - grnd()*((v1+v2)**3 - abs(v1-v2)**3))**(2.0d0/3.0d0))/vr
    vb = sqrt(max(0.0d0, (2*v1*v2/vr)**2 - va**2))
    
    dx = vr/2 - va
    vcs = dx**2 + vb**2
    
    vcs1 = sqrt(vcs)
    ran = sqrt(grnd())
    ran1 = sqrt(1 - ran**2)
    ran2 = grnd()*2*pi
    tau = ran1*cos(ran2)
    ex = ran*dx/vcs1 - tau*vb/vcs1
    ey = ran*vb/vcs1 + tau*dx/vcs1
    ez = ran1*sin(ran2)
    
    if (BounceIf(dx, vb, 0.0d0, ex, ey, ez, i2, j2)) then
      e = ecoef(dx, vb, 0.0d0, ex, ey, ez, i2, j2)
      vcs = (1+e)*vcs1*ran
      
      if (v1 < v2) then
        tmpi = i1
        tmpi2 = j1
        i1 = tmpi2
        j1 = tmpi
        tmpi = i2
        tmpi2 = j2
        i2 = tmpi2
        j2 = tmpi
        tmp = v1
        tmp2 = v2
        v1 = tmp2
        v2 = tmp
      end if
      
      vcs1 = vcs*j2/(i2 + j2)
      vcs2 = vcs*i2/(i2 + j2)
      vcs1 = sqrt(v1**2 + vcs1**2 - 2*vcs1*v1*ex)
      vcs2 = sqrt(v2**2 + vcs2**2 + 2*vcs2*(va*ex - vb*ey))
      vs%d(i1) = vcs1
      if (vcs1 > vmaxvec%d(i2)) then
        vmaxvec%d(i2) = vcs1
        uioi%d(:) = ui%d(i2,:)
        vioi%d(:) = vi%d(i2,:)
        	
        call OneUpdate(i2, u, v, ui, vi, partofsize, vmaxvec)
        	
        uioi%d(:) = ui%d(i2,:) - uioi%d(:)
        vioi%d(:) = vi%d(i2,:) - vioi%d(:)
        	
        call treeupdate2(utree, ui, i2, treestart, uioi)
        call treeupdate2(vtree, vi, i2, treestart, vioi)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
      end if
      vs%d(j1) = vcs2
      if (vcs2 > vmaxvec%d(j2)) then
        vmaxvec%d(j2) = vcs2
        uioi%d(:) = ui%d(j2,:)
        vioi%d(:) = vi%d(j2,:)
        	
        call OneUpdate(j2, u, v, ui, vi, partofsize, vmaxvec)
        	
        uioi%d(:) = ui%d(j2,:) - uioi%d(:)
        vioi%d(:) = vi%d(j2,:) - vioi%d(:)
        	
        call treeupdate2(utree, ui, j2, treestart, uioi)
        call treeupdate2(vtree, vi, j2, treestart, vioi)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
      end if
    else
      !Needed for right order remove in partofnum
      if (i0 > j0) then
        tmpi = i0
        tmpi2 = j0
        i0 = tmpi2
        j0 = tmpi
        tmpi = i1
        tmpi2 = j1
        i1 = tmpi2
        j1 = tmpi
        tmpi = i2
        tmpi2 = j2
        i2 = tmpi2
        j2 = tmpi
        tmp = v1
        tmp2 = v2
        v1 = tmp2
        v2 = tmp
      end if
      partofnum(j2)%d(j0) = partofnum(j2)%d(partofsize%d(j2))
      partofsize%d(j2) = partofsize%d(j2) - 1
      partofnum(i2)%d(i0) = partofnum(i2)%d(partofsize%d(i2))
      partofsize%d(i2) = partofsize%d(i2) - 1
      k2 = i2 + j2
      if (k2 > maxmass) then
        Allocate(tmpintvec(maxmass))
        tmpintvec(:) = partofnum(:)
        Deallocate(partofnum)
        Allocate(partofnum(4*maxmass))
        partofnum(1:maxmass) = tmpintvec(1:maxmass)
        do j = 1+maxmass,4*maxmass
          call partofnum(j)%init(partmin)
        end do
        Deallocate(tmpintvec)
        Allocate(tmpint(maxmass))
        tmpint(:) = partofsize%d(:)
        call partofsize%deinit()
        call partofsize%init(4*maxmass)
        partofsize%d(1:maxmass) = tmpint(1:maxmass)
        Deallocate(tmpint)
        call tmpvec%init(maxmass)
        tmpvec%d(:) = vmaxvec%d(:)
        call vmaxvec%deinit()
        call vmaxvec%init(4*maxmass)
        vmaxvec%d(1:maxmass) = tmpvec%d(1:maxmass)
        call tmpvec%deinit()
        maxmass = maxmass*4
        
        call u%deinit()
        call v%deinit()
        call ui%deinit()
        call vi%deinit()
        call FullApproximate(u, v, ui, vi, partofsize, vmaxvec)
        treestart = treestart + maxmass/4
        utree = treebuild(ui, treestart)
        vtree = treebuild(vi, treestart)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
        
      end if
      if (partofsize%d(k2) >= partofnum(k2)%n) then
        Allocate(tmpint(partofsize%d(k2)))
        tmpint(:) = partofnum(k2)%d(:)
        call partofnum(k2)%deinit()
        call partofnum(k2)%init(2*partofsize%d(k2))
        partofnum(k2)%d(1:partofsize%d(k2)) = tmpint(1:partofsize%d(k2))
        Deallocate(tmpint)
      end if
      partofsize%d(k2) = partofsize%d(k2) + 1
      partofnum(k2)%d(partofsize%d(k2)) = i1

      vmass%d(i1) = k2
      vmass%d(j1) = 0
      vcs = sqrt((i2*v1/k2)**2 + (j2*v2/k2)**2 + dble(i2)*j2*vr/k2*va/k2)
      vs%d(i1) = vcs
      vs%d(j1) = 0
      if (vcs > vmaxvec%d(k2)) then
        vmaxvec%d(k2) = vcs
      end if
      
      uioi%d(:) = ui%d(i2,:)
      vioi%d(:) = vi%d(i2,:)
      
      call OneUpdate(i2, u, v, ui, vi, partofsize, vmaxvec)
      
      uioi%d(:) = ui%d(i2,:) - uioi%d(:)
      vioi%d(:) = vi%d(i2,:) - vioi%d(:)
      
      call treeupdate2(utree, ui, i2, treestart, uioi)
      call treeupdate2(vtree, vi, i2, treestart, vioi)
      
      uioj%d(:) = ui%d(j2,:)
      vioj%d(:) = vi%d(j2,:)
      
      call OneUpdate(j2, u, v, ui, vi, partofsize, vmaxvec)
      
      uioj%d(:) = ui%d(j2,:) - uioj%d(:)
      vioj%d(:) = vi%d(j2,:) - vioj%d(:)
      
      call treeupdate2(utree, ui, j2, treestart, uioj)
      call treeupdate2(vtree, vi, j2, treestart, vioj)
      
      uiok%d(:) = ui%d(k2,:)
      viok%d(:) = vi%d(k2,:)
      
      call OneUpdate(k2, u, v, ui, vi, partofsize, vmaxvec)
      
      uiok%d(:) = ui%d(k2,:) - uiok%d(:)
      viok%d(:) = vi%d(k2,:) - viok%d(:)
      
      call treeupdate2(utree, ui, k2, treestart, uiok)
      call treeupdate2(vtree, vi, k2, treestart, viok)
      
      usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
      vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
      totalvec%d(:) = usum%d(:)*vsum%d(:)
      totalrate = sum(totalvec%d)
      
      np = np - 1
      
      !Number of particles doubling
      if (np <= np0/2) then
        np = np*2
        i = 1
        do j = 1, maxmass
          j1 = partofnum(j)%n
          if (2*partofsize%d(j) > j1) then
            Allocate(tmpint(partofsize%d(j)))
            tmpint(1:partofsize%d(j)) = partofnum(j)%d(1:partofsize%d(j))
            call partofnum(j)%deinit()
            call partofnum(j)%init(2*j1)
            partofnum(j)%d(1:partofsize%d(j)) = tmpint(1:partofsize%d(j))
            Deallocate(tmpint)
          end if
          do j1 = 1 + partofsize%d(j), 2*partofsize%d(j)
            do while(vmass%d(i) .ne. 0)
              i = i + 1
            end do
            vmass%d(i) = j
            j2 = partofnum(j)%d(j1-partofsize%d(j))
            vs%d(i) = vs%d(j2)
            partofnum(j)%d(j1) = i
          end do
        end do
        partofsize%d(:) = partofsize%d(:)*2
        ui = ui*2.0d0
        vi = vi*2.0d0
        utree = treebuild(ui, treestart)
        vtree = treebuild(vi, treestart)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
        n0 = n0/2
      end if
    end if
    
    !Thermostat
    if ((mod(step,Nht) == 0) .and. (gbig > 0)) then
      tau = gbig*(curtime-lasttime)
      do j2 = 1, maxmass
        if (Tmolgas > 0) then
          tau2 = sqrt(3*Tmolgas/j2*(1 - exp(-tau*(j2**(gsmall+1))/Tmolgas)))
          tauexp = exp(-tau*(j2**(gsmall+1))/Tmolgas/2)
        else
          tau2 = sqrt(3*tau*j2**(gsmall))
          tauexp = 1.0d0
        end if
        do j1 = 1, partofsize%d(j2)
          j = partofnum(j2)%d(j1)
          call randiradv(ex, ey, ez)
          !Rewrite
          !vx%d(j) = vx%d(j)*tauexp+ex*tau2
          !vy%d(j) = vy%d(j)*tauexp+ey*tau2
          !vz%d(j) = vz%d(j)*tauexp+ez*tau2
        end do
      end do
      lasttime = curtime
    
      if (mod(step,Nh) .ne. 0) then
        vmaxvec%d = 0
        do j = 1, np0
          j1 = vmass%d(j)
          if (j1 > 0) then
            vmaxvec%d(j1) = max(vmaxvec%d(j1), vs%d(j))
          end if
        end do
      
        call FullUpdate(u, v, ui, vi, partofsize, vmaxvec)
        utree = treebuild(ui, treestart)
        vtree = treebuild(vi, treestart)
        usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
        vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
        totalvec%d(:) = usum%d(:)*vsum%d(:)
        totalrate = sum(totalvec%d)
      end if
    end if
    
    if ((mod(step,Nh) == 0) .or. (curtime >= maxtime)) then
    
      vmaxvec%d = 0
      tempsum = 0
      call temps%deinit()
      call temps%init(maxmass)
      do j = 1, np0
        j1 = vmass%d(j)
        if (j1 > 0) then
          vmaxvec%d(j1) = max(vmaxvec%d(j1), vs%d(j))
          temps%d(j1) = temps%d(j1) + vs%d(j)**2*j1
          tempsum = tempsum + vs%d(j)**2*j1
        end if
      end do
      ones = eveci(maxmass)
      temps%d(1:maxmass) = temps%d(1:maxmass)/max(ones%d(1:maxmass),partofsize%d(1:maxmass))/3
      if (verbose > 1) then
        print *, 'TOTAL ENERGY', tempsum
      end if
      tavg = tempsum/np/3
      
      call FullUpdate(u, v, ui, vi, partofsize, vmaxvec)
      utree = treebuild(ui, treestart)
      vtree = treebuild(vi, treestart)
      usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
      vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
      totalvec%d(:) = usum%d(:)*vsum%d(:)
      totalrate = sum(totalvec%d)

      if (verbose > 1) then
!         write(2, *) curtime, temps%d(1), temps%d(3), temps%d(5), temps%d(8), tavg
!         write(3, *) curtime, dble(partofsize%d(1))/np0, (n0*np)/np0
      
        print *, step, curtime, np, maxmass
        print *, 'Temperature', temps%d(1), temps%d(2), tavg
        print *, 'Time passed', dsecnd()-cursec
        
        timecur = curtime
        timeold = timecur
        
      end if
    end if
  end do
!   close(2)
!   close(3)
  n0 = (n0*np)/np0
  if (verbose > 0) then
    if (verbose == 1) then
      print *, 'Time passed', dsecnd()-cursec
    end if
    print *, 'Collisions', step
    print *, 'Number density', n0
    print *, 'DONE'
  end if
  
  Deallocate(sizestart)
  do j = 1, size(partofnum)
    call partofnum(j)%deinit()
  end do
  Deallocate(partofnum)
end

subroutine DSMCup(n0, partofsize, temps, maxtime, gbig, gsmall, Tmolgas, verbose)
  Double precision, intent(inout) :: n0 !Total number density; IN and OUT
  Type(IntVec) :: partofsize !Numbers of particles; IN and OUT
  Type(Vector) :: temps !Temperatures of each size; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Double precision, intent(in) :: gbig, gsmall, tmolgas !Thermostat parameters
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every Nh steps
  
  Integer(4) np, np0 !Current and initial total number of particles
  Double precision e !Restitution coefficient
  Type(Vector) vx, vy, vz !Particle speeds
  Type(Vector) vs !Particle squared speed
  Type(Vector) vmaxvec !Maximum speed for each size
  Integer(4), allocatable :: sizestart(:) !Start index of size-i particles in the initial array
  Integer(4), allocatable :: vmass(:) !Integer particle masses
  
  Double precision curtime !Current system time
  Double precision ex, ey, ez !Collision direction components
  Double precision vcx, vcy, vcz, vcs, vcs1, vcs2 !Temporaries for speed values
  Double precision cursec !Current time
  Double precision dsecnd !LAPACK procedure for time measurement
  Integer(4) maxmass !Maximum cluster mass in the system
  
  Integer(8) Nh !Update output every Nh steps
  Integer(8) step !Current step
  Double precision tempsum !Total kinetic energy (sum of weighted temperatures)
  Integer(4) i, j, i1, j1, tmpi, i2, j2 !Indices
  
  Double precision vmin !Minimum accepted speed
  Double precision cijmax !Collision kernel upper bound
  
  !Temporaries
  Double precision tmp, tmp2
  
  maxmass = partofsize%n
  np = sum(partofsize%d)
  np0 = np
  
  !Update time
  Nh = np!/20!5
  
  Allocate(sizestart(maxmass+1))
  Allocate(vmass(np))
  sizestart(1) = 1
  i = 1
  j1 = 1
  do while (i <= maxmass)
    tmpi = partofsize%d(i)
    do j = j1, j1+tmpi-1
      vmass(j) = i
    end do
    j1 = j1 + tmpi
    i = i + 1
    sizestart(i) = sizestart(i-1) + tmpi
  end do

  call vx%init(np)
  call vy%init(np)
  call vz%init(np)
  call vs%init(np+1)
  
  do i = 1, np
    vx%d(i) = ggrnd()
    vy%d(i) = ggrnd()
    vz%d(i) = ggrnd()
  end do
  
  do j = 1, maxmass
    ex = 0.0d0
    ey = 0.0d0
    ez = 0.0d0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      ex = ex + vx%d(i)
      ey = ey + vy%d(i)
      ez = ez + vz%d(i)
    end do
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vx%d(i) = vx%d(i) - ex/partofsize%d(j)
      vy%d(i) = vy%d(i) - ey/partofsize%d(j)
      vz%d(i) = vz%d(i) - ez/partofsize%d(j)
    end do
    tmp = 0.0d0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      tmp = tmp + vx%d(i)**2 + vy%d(i)**2 + vz%d(i)**2
    end do
    tmp = tmp*j
    if (partofsize%d(j) > 1) then
      tmp = 3*temps%d(j)/tmp*partofsize%d(j)
    elseif (partofsize%d(j) == 1) then
      call randir(vx%d(sizestart(j)), vy%d(sizestart(j)), vz%d(sizestart(j)))
      tmp = 3*temps%d(j)/j
    end if
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vx%d(i) = vx%d(i)*sqrt(tmp)
      vy%d(i) = vy%d(i)*sqrt(tmp)
      vz%d(i) = vz%d(i)*sqrt(tmp)
    end do
  end do
  
  do i = 1, np
    vx%d(i) = vx%d(i)/sqrt(temps%d(vmass(i)))
    vy%d(i) = vy%d(i)/sqrt(temps%d(vmass(i)))
    vz%d(i) = vz%d(i)/sqrt(temps%d(vmass(i)))
  end do
  do i = 1, np
    vs%d(i) = vx%d(i)**2+vy%d(i)**2+vz%d(i)**2
  end do
  tempsum = 0
  do j = 1, maxmass
    temps%d(j) = 0.0d0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      temps%d(j) = temps%d(j) + vs%d(i)
    end do
    tempsum = tempsum + temps%d(j)*j
    temps%d(j) = (temps%d(j)*j)/max(1,partofsize%d(j))/3
  end do
  if (verbose > 0) then
    print *, 'Monomer Temperature', temps%d(1)
    print *, 'TOTAL ENERGY START', tempsum
  end if
  
  call vmaxvec%init(maxmass)
  do j = 1, maxmass
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vmaxvec%d(j) = max(vmaxvec%d(j), vs%d(i))
    end do
    vmaxvec%d(j) = sqrt(vmaxvec%d(j))
  end do
  
  cijmax = 0
  do i = 1, np
    cijmax = max(cijmax, 2*(vmass(i)**(1.0d0/3.0d0) + maxmass**(1.0d0/3.0d0))**2*sqrt(vs%d(i)))
  end do
  
  curtime = 0.0d0
  step = 0
!   open(2, file='DSMCupTime')
  cursec = dsecnd()
  do while (curtime < maxtime)
    step = step + 1
    i = 0
    do while (1 > 0)
      curtime = curtime - log(grnd())/n0*np0/(np*(np-1)*cijmax)/pi
      
      i1 = np+1
      j1 = np+1
      do while (i1 > np)
        i1 = floor(grnd()*np)+1
      end do
      do while ((j1 > np) .or. (j1 == i1))
        j1 = floor(grnd()*np)+1
      end do
      i2 = vmass(i1)
      j2 = vmass(j1)
      
      call randir(ex, ey, ez)
      vmin = grnd()*cijmax
      if (i1 .ne. j1) then
        vcx = ex*(vx%d(i1)-vx%d(j1))
        vcy = ey*(vy%d(i1)-vy%d(j1))
        vcz = ez*(vz%d(i1)-vz%d(j1))
        vcs = vcx+vcy+vcz
        if (abs(vcs)*(i2**(1.0d0/3.0d0) + j2**(1.0d0/3.0d0))**2 > vmin) then
          exit
        end if
        if ((sqrt(vs%d(i1))+sqrt(vs%d(j1)))*(i2**(1.0d0/3.0d0) + j2**(1.0d0/3.0d0))**2 > cijmax) then
          print *, i2, j2, maxmass
          stop 0
        end if
      end if
    end do
    e = 0.99d0
    if ((e > 1) .or. (e < 0)) then
      print *, 'E ERROR', e
    end if
    vcs = (1+e)*vcs
    vcs1 = vcs*vmass(j1)/dble(vmass(i1) + vmass(j1))
    vcs2 = vcs*vmass(i1)/dble(vmass(i1) + vmass(j1))
    vx%d(i1) = vx%d(i1) - vcs1*ex
    vy%d(i1) = vy%d(i1) - vcs1*ey
    vz%d(i1) = vz%d(i1) - vcs1*ez
    vx%d(j1) = vx%d(j1) + vcs2*ex
    vy%d(j1) = vy%d(j1) + vcs2*ey
    vz%d(j1) = vz%d(j1) + vcs2*ez
    
    vs%d(i1) = vx%d(i1)**2+vy%d(i1)**2+vz%d(i1)**2
    vs%d(j1) = vx%d(j1)**2+vy%d(j1)**2+vz%d(j1)**2
    
    cijmax = max(cijmax, 2*(vmass(i1)**(1.0d0/3.0d0) + maxmass**(1.0d0/3.0d0))**2*sqrt(vs%d(i1)))
    cijmax = max(cijmax, 2*(vmass(j1)**(1.0d0/3.0d0) + maxmass**(1.0d0/3.0d0))**2*sqrt(vs%d(j1)))
    
    if ((mod(step,Nh) == 0) .or. (curtime >= maxtime)) then
    
      tmp = sum(vmass)
      tmp2 = sum(vx%d * vmass)
      ex = tmp2/tmp
      tmp2 = sum(vy%d * vmass)
      ey = tmp2/tmp
      tmp2 = sum(vz%d * vmass)
      ez = tmp2/tmp
      do i = 1, np
        vx%d(i) = vx%d(i) - ex
        vy%d(i) = vy%d(i) - ey
        vz%d(i) = vz%d(i) - ez
      end do
      
      vmaxvec%d = 0
      tempsum = 0
      do j = 1, maxmass
        vmaxvec%d(j) = 0
        temps%d(j) = 0
        do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
          vs%d(i) = vx%d(i)**2 + vy%d(i)**2 + vz%d(i)**2
          temps%d(j) = temps%d(j) + vs%d(i)
          vmaxvec%d(j) = max(vmaxvec%d(j), vs%d(i))
        end do
        vmaxvec%d(j) = sqrt(vmaxvec%d(j))
        tempsum = tempsum + temps%d(j)*j
        temps%d(j) = (temps%d(j)*j)/partofsize%d(j)/3
      end do
      if (verbose > 1) then
        print *, 'TOTAL ENERGY', tempsum
      end if
      
      cijmax = 0
      do i = 1, np
        cijmax = max(cijmax, 2*(vmass(i)**(1.0d0/3.0d0) + maxmass**(1.0d0/3.0d0))**2*sqrt(vs%d(i)))
      end do
    
      if (verbose > 1) then
!         write(2, *) curtime, temps%d(1), temps%d(2), temps%d(5), temps%d(10), temps%d(20), temps%d(50)
        print *, step, curtime
        print *, 'Temperature', temps%d(1), temps%d(2)
        print *, 'Time passed', dsecnd()-cursec
      end if
    end if
  end do
  tempsum = 0
  do j = 1, maxmass
    vmaxvec%d(j) = 0
    temps%d(j) = 0
    do i = sizestart(j), partofsize%d(j)+sizestart(j)-1
      vs%d(i) = vx%d(i)**2 + vy%d(i)**2 + vz%d(i)**2
      temps%d(j) = temps%d(j) + vs%d(i)
      vmaxvec%d(j) = max(vmaxvec%d(j), vs%d(i))
    end do
    vmaxvec%d(j) = sqrt(vmaxvec%d(j))
    tempsum = tempsum + temps%d(j)*j
    temps%d(j) = (temps%d(j)*j)/partofsize%d(j)/3
  end do
  print *, 'Time passed', dsecnd()-cursec
  print *, 'Collisions', step
  print *, 'Tavg', tempsum/np/3
  print *, 'Number density', n0
  print *, 'DONE'
  Deallocate(sizestart)
  Deallocate(vmass)
end

!Initial DSMC kernel approximation
subroutine FullApproximate(u, v, ui, vi, nvec, vvec)
  Type(Mtrx), intent(out) :: u, v, ui, vi
  Type(IntVec), intent(in) :: nvec
  Type(Vector), intent(in) :: vvec

  Integer(4) n, r
  
  n = nvec%n
  call DSMCKernel(n, u, v)
  r = u%m
  
  call ui%init(n, r)
  call vi%init(n, r)
  call FullUpdate(u, v, ui, vi, nvec, vvec)
end

!Update of DSMC kernel using maximum speed values
subroutine FullUpdate(u, v, ui, vi, nvec, vvec)
  Type(Mtrx), intent(in) :: u, v
  Type(Mtrx) :: ui, vi
  Type(IntVec), intent(in) :: nvec
  Type(Vector), intent(in) :: vvec
  
  Integer(4) n, i, j, r
  
  !ui * vi = C_{ij} = sigma_{ij} * (vmax(i) + vmax(j)) n_i n_j  ~ (2 sigma_{ij}) * vmax(i) n_i n_j
  n = u%n
  r = u%m
  do i = 1, r
    do j = 1, n
      ui%d(j,i) = u%d(j,i)*Uspeed(vvec%d(j), j, i)*nvec%d(j)
      vi%d(j,i) = v%d(j,i)*Vspeed(vvec%d(j), j, i)*nvec%d(j)
    end do
  end do
  
end

!Update i-th low-rank factors in DSMC kernel
subroutine OneUpdate(i, u, v, ui, vi, nvec, vvec)
  Integer(4), intent(in) :: i
  Type(Mtrx), intent(in) :: u, v
  Type(Mtrx) :: ui, vi
  Type(IntVec), intent(in) :: nvec
  Type(Vector), intent(in) :: vvec
  
  Integer(4) j

  do j = 1, u%m
    ui%d(i,j) = u%d(i,j)*Uspeed(vvec%d(i), i, j)*nvec%d(i)
    vi%d(i,j) = v%d(i,j)*Vspeed(vvec%d(i), i, j)*nvec%d(i)
  end do
end

!Pregenerate random numbers?
subroutine MonteTemp(nv0, nvi, temps, maxtime, lambda, apar, verbose)
  Double precision, intent(inout) :: nv0 !Total number density; IN and OUT
  Type(IntVec) :: nvi !Numbers of particles; IN and OUT
  Type(Vector) :: temps !Temperatures of each size; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Double precision, intent(in) :: lambda, apar !Kernel parameters
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  
  Integer(4) np, maxnp !Current and initial total number of particles
  Integer(4) maxmass !Maximum cluster mass in the system
  Integer(4) treestart !Start index of segment tree data
  Integer(4) r, curr !Rank, iteration over ranks
  Integer(4) posit !True if low-rank factors are positive
  Integer(8) step !Current step
  Integer(8) everys !How often write to file
  Integer(8) everys2 !How often apply thermostat
  Double precision nv1 !Number density per particle
  Double precision tavg !Current average temperature
  Double precision prevcurtime, curtime !Current and previous thermostat system times
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  Double precision randi !Random number
  Double precision totalrate !Total collision rate
  Type(Mtrx) u, v, ui, vi !Low-rank factors
  Type(Mtrx) utree, vtree !Segment trees
  Type(Vector) usum, vsum, totalvec !Vectors, containing sums of collision rates
  Double precision gsmall, Tmolgas !Molecular gas parameters
  Double precision tauexp, tauadd !Thermostat changes
  
  Integer(4) j, i1, j1, k1 !Indices
  
  !Temporaries
  Double precision tau, tmp1, tmp2
  Type(Vector) tmpvec, uioi, uioj, uiok, vioi, vioj, viok
  Type(IntVec) nvc
  
  maxmass = max(16, ibset(0, 2*((bit_size(nvi%n-1) - leadz(nvi%n-1) + 1)/2)))
  call nvc%copy(nvi)
  call nvi%deinit()
  call nvi%init(maxmass)
  nvi%d(1:nvc%n) = nvc%d(1:nvc%n)
  call nvc%deinit()
  call tmpvec%copy(temps)
  call temps%deinit()
  call temps%init(maxmass)
  temps%d(1:tmpvec%n) = tmpvec%d(1:tmpvec%n)
  call tmpvec%deinit()
  
  !Molecular gas effect
  gsmall = -1.0d0/3.0d0 !(-4/3)
  Tmolgas = 0.0d0!1.0d0

  treestart = 0
  !Calculate geometric series
  i1 = 1
  do while (i1*4 < maxmass)
    treestart = treestart + i1
    i1 = i1*4
  end do
  
  np = sum(nvi%d)
  maxnp = np
  everys = np!np/5
  everys2 = np!np/100
  
  nv1 = nv0/np

  call TempApproximate(u, v, ui, vi, nvi, temps, posit)
  r = u%m
  call usum%init(r)
  call vsum%init(r)
  call uioi%init(r)
  call uioj%init(r)
  call uiok%init(r)
  call vioi%init(r)
  call vioj%init(r)
  call viok%init(r)
  call totalvec%init(r)
  utree = treebuild(ui, treestart)
  vtree = treebuild(vi, treestart)
  usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
  vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
  totalvec%d(:) = usum%d(:)*vsum%d(:)
  totalrate = sum(totalvec%d)
  
!   open(2, file='TimeTemp')
  curtime = 0.0d0
  step = 0
  prevcurtime = curtime
  dstart = dsecnd()
  do while (curtime < maxtime)
    step = step + 1

    i1 = 0
    do while (1 > 0)
    
      tau = 1.0d0/totalrate/nv1
      curtime = curtime + tau
      
      if (posit == 1) then
        randi = grnd()*totalrate
        do curr = 1, r-1
          randi = randi - totalvec%d(curr)
          if (randi <= 0) then
            exit
          end if
        end do
        
        randi = grnd()*usum%d(curr)
        i1 = treefindr(utree, ui, randi, curr, treestart)
        randi = grnd()*vsum%d(curr)
        j1 = treefindr(vtree, vi, randi, curr, treestart)
      else
        randi = grnd()*totalrate
        i1 = treefind(utree, ui, vsum, randi, treestart)
        tmpvec = tovec(ui%subarray(i1,r,i1,1))
        randi = grnd()*(tmpvec*vsum)
        j1 = treefind(vtree, vi, tmpvec, randi, treestart)
        call tmpvec%deinit()
      end if
      
      if (grnd() < sucrate(i1, j1, temps%d(i1), temps%d(j1))) then
        if (i1 .ne. j1) then
          exit
        else
          if (1.0d0 < grnd() * nvi%d(i1)) then
            exit
          end if
        end if
      end if
      
    end do

    if ((temps%d(i1) + temps%d(j1) > 0) .and. (nvi%d(i1) .ne. 0) .and. (nvi%d(j1) .ne. 0)) then
    
      if (grnd() < BounceProb(lambda, apar, i1, j1, temps%d(i1), temps%d(j1))) then
      
        tmp1 = BounceMinus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), nvi%d(i1))
        tmp2 = BounceMinus(lambda, apar, j1, i1, temps%d(j1), temps%d(i1), nvi%d(j1))
        if (i1 == j1) then
          tmp2 = 2*tmp1 - temps%d(j1)
          tmp1 = tmp2
        end if
        temps%d(i1) = tmp1
        temps%d(j1) = tmp2
        
        uioi%d(1:r) = ui%d(i1,1:r)
        uioj%d(1:r) = ui%d(j1,1:r)
        vioi%d(1:r) = vi%d(i1,1:r)
        vioj%d(1:r) = vi%d(j1,1:r)
        
        call TempUpdate(i1, u, v, ui, vi, nvi, temps)
        call TempUpdate(j1, u, v, ui, vi, nvi, temps)
        
        uioi%d(1:r) = ui%d(i1,1:r) - uioi%d(1:r)
        uioj%d(1:r) = ui%d(j1,1:r) - uioj%d(1:r)
        vioi%d(1:r) = vi%d(i1,1:r) - vioi%d(1:r)
        vioj%d(1:r) = vi%d(j1,1:r) - vioj%d(1:r)
        if (i1 == j1) then
          uioj%d = 0
          vioj%d = 0
        end if
      
        call treeupdate2(utree, ui, i1, treestart, uioi)
        call treeupdate2(utree, ui, j1, treestart, uioj)
        call treeupdate2(vtree, vi, i1, treestart, vioi)
        call treeupdate2(vtree, vi, j1, treestart, vioj)
        
      else
      
        k1 = i1 + j1
        
        if (k1 > maxmass) then
        
          call tmpvec%copy(temps)
          call temps%deinit()
          call temps%init(maxmass*4)
          temps%d(1:maxmass) = tmpvec%d(1:maxmass)
          temps%d(maxmass+1:maxmass*4) = 1
          call tmpvec%deinit()
          Allocate(nvc%d(maxmass))
          nvc%d(1:maxmass) = nvi%d(1:maxmass)
          call nvi%deinit()
          call nvi%init(maxmass*4)
          nvi%d(1:maxmass) = nvc%d(1:maxmass)
          Deallocate(nvc%d)
          
          call u%deinit()
          call v%deinit()
          call ui%deinit()
          call vi%deinit()
          call TempApproximate(u, v, ui, vi, nvi, temps, posit)
          call utree%deinit()
          call vtree%deinit()
          treestart = treestart + maxmass/4
          utree = treebuild(ui, treestart)
          vtree = treebuild(vi, treestart)
          maxmass = maxmass*4
          
        end if
        
        tmp1 = AggloMinus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), nvi%d(i1))
        tmp2 = AggloMinus(lambda, apar, j1, i1, temps%d(j1), temps%d(i1), nvi%d(j1))
        if (i1 == j1) then
          tmp2 = max(0.0d0, 2*tmp1 - temps%d(j1))
          tmp1 = tmp2
        end if
        temps%d(k1) = AggloPlus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), temps%d(k1), nvi%d(k1))
        temps%d(i1) = tmp1
        temps%d(j1) = tmp2
        
        np = np-1
        nvi%d(i1) = nvi%d(i1)-1
        nvi%d(j1) = nvi%d(j1)-1
        nvi%d(k1) = nvi%d(k1)+1
        
        uioi%d(1:r) = ui%d(i1,1:r)
        uioj%d(1:r) = ui%d(j1,1:r)
        uiok%d(1:r) = ui%d(k1,1:r)
        vioi%d(1:r) = vi%d(i1,1:r)
        vioj%d(1:r) = vi%d(j1,1:r)
        viok%d(1:r) = vi%d(k1,1:r)
        
        call TempUpdate(k1, u, v, ui, vi, nvi, temps)
        call TempUpdate(i1, u, v, ui, vi, nvi, temps)
        call TempUpdate(j1, u, v, ui, vi, nvi, temps)
        
        uioi%d(1:r) = ui%d(i1,1:r) - uioi%d(1:r)
        uioj%d(1:r) = ui%d(j1,1:r) - uioj%d(1:r)
        uiok%d(1:r) = ui%d(k1,1:r) - uiok%d(1:r)
        vioi%d(1:r) = vi%d(i1,1:r) - vioi%d(1:r)
        vioj%d(1:r) = vi%d(j1,1:r) - vioj%d(1:r)
        viok%d(1:r) = vi%d(k1,1:r) - viok%d(1:r)
        
        if (i1 == j1) then
          uioj%d = 0
          vioj%d = 0
        end if
        
        call treeupdate2(utree, ui, k1, treestart, uiok)
        call treeupdate2(utree, ui, i1, treestart, uioi)
        call treeupdate2(utree, ui, j1, treestart, uioj)
        call treeupdate2(vtree, vi, k1, treestart, viok)
        call treeupdate2(vtree, vi, i1, treestart, vioi)
        call treeupdate2(vtree, vi, j1, treestart, vioj)
        
      end if
      
      usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
      vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
      totalvec%d(:) = usum%d(:)*vsum%d(:)
      totalrate = sum(totalvec%d)
      
    end if

    if (np <= maxnp/2) then
      np = np*2
      nvi%d(:) = nvi%d(:)*2
      ui = ui*2.0d0
      vi = vi*2.0d0
      utree = treebuild(ui, treestart)
      vtree = treebuild(vi, treestart)
      usum%d(:) = utree%d(1,:) + utree%d(2,:) + utree%d(3,:) + utree%d(4,:)
      vsum%d(:) = vtree%d(1,:) + vtree%d(2,:) + vtree%d(3,:) + vtree%d(4,:)
      totalvec%d(:) = usum%d(:)*vsum%d(:)
      totalrate = sum(totalvec%d)
      
      nv1 = nv1/2
    end if
    
    if ((mod(step,everys2) == 0) .or. (curtime >= maxtime)) then
    
      if (Tmolgas > 0) then
        tau = curtime - prevcurtime
        do j = 1, maxmass
          tauexp = exp(-tau/j**(1.0d0/3.0d0)/Tmolgas)
          tauadd = Tmolgas*(1 - tauexp)
          temps%d(j) = temps%d(j)*tauexp + tauadd
        end do
        prevcurtime = curtime
      end if
    
      tavg = 0.0d0
      do j = 1, maxmass
        tavg = tavg + temps%d(j)*nvi%d(j)
      end do
      tavg = tavg/np
      
      if (verbose > 1) then
!         write(2, *) curtime, temps%d(1), temps%d(3), temps%d(5), temps%d(10), temps%d(20), tavg,  nv1*np
      end if
    end if
    if ((verbose > 1) .and. ((mod(step,everys) == 0) .or. (curtime >= maxtime))) then
      print *, step, curtime, np, maxmass
      print *, 'Temperatures', temps%d(1), temps%d(2), tavg
      print *, 'Time passed', dsecnd() - dstart
    end if
  end do
!   close(2)
  nv0 = nv1*np
  
  if (verbose > 0) then
    if (verbose == 1) then
      print *, 'Time passed', dsecnd() - dstart
    end if
    print *, 'Collisions', step
  
    print *, 'Number density', nv0
    tau = 0
    do j = 1, maxmass
      tau = tau + nv1*dble(nvi%d(j))*dble(j)*dble(j)
    end do
    print *, 'Second moment, M2', tau
    print *, 'Tavg', tavg
    print *, 'Maxmass', maxmass
    print *, 'DONE'
  end if
end

!maxmass should be a power of 4
!fastest trees.
subroutine MonteSimple(nv0, nvi, maxtime, verbose)
  Double precision, intent(inout) :: nv0 !Total number density; IN and OUT
  Type(IntVec) :: nvi !Numbers of particles; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  
  Integer(4) np, maxnp !Current and initial total number of particles
  Integer(4) maxmass !Maximum cluster mass in the system
  Integer(4) treestart !Start index of segment tree data
  
  Integer(8) step !Current step
  Integer(8) everys !How often write to file
  Double precision nv1 !Number density per particle
  Double precision curtime !Current system time
  Double precision randi !Random number
  Double precision totalrate !Total collision rate
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  
  Type(Vector) u, v, ui, vi !Low-rank factors
  Type(Vector) utree, vtree !Segment trees
  Double precision usum, vsum !Sums of collision rates
  
  Integer(4) j, i1, j1, k1 !Indices
  
  !Temporaries
  Double precision tau
  Type(IntVec) nvc
  
  maxmass = max(16, ibset(0, 2*((bit_size(nvi%n-1) - leadz(nvi%n-1) + 1)/2)))
  call nvc%copy(nvi)
  call nvi%deinit()
  call nvi%init(maxmass)
  nvi%d(1:nvc%n) = nvc%d(1:nvc%n)
  call nvc%deinit()

  treestart = 0
  !Calculate geometric series
  i1 = 1
  do while (i1*4 < maxmass)
    treestart = treestart + i1
    i1 = i1*4
  end do
  
  np = sum(nvi%d)
  maxnp = np
  everys = np!/5
  
  nv1 = nv0/np

  call u%init(maxmass)
  call v%init(maxmass)
  call ui%init(maxmass)
  call vi%init(maxmass)
  do j = 1, maxmass
    u%d(j) = u_vec(j)
    v%d(j) = v_vec(j)
  end do
  do j = 1, maxmass
    ui%d(j) = nvi%d(j)*u%d(j)
    vi%d(j) = nvi%d(j)*v%d(j)
  end do
  utree = treebuildv(ui, treestart)
  usum = utree%d(1) + utree%d(2) + utree%d(3) + utree%d(4)
  vtree = treebuildv(vi, treestart)
  vsum = vtree%d(1) + vtree%d(2) + vtree%d(3) + vtree%d(4)
  totalrate = usum*vsum
  
!   open(2, file='TimeLowrank')
  curtime = 0.0d0
  step = 0
  dstart = dsecnd()
  do while (curtime < maxtime)
    step = step + 1

    i1 = 0
    do while (1 > 0)
    
      tau = 1.0d0/totalrate/nv1
      curtime = curtime + tau
      
      randi = grnd()*usum
      i1 = treefindv(utree, ui, randi, treestart)
      randi = grnd()*vsum
      j1 = treefindv(vtree, vi, randi, treestart)
      
      if (grnd() < sucrates(i1, j1)) then
        if (i1 .ne. j1) then
          exit
        else
          if (1.0d0 < grnd() * nvi%d(i1)) then
            exit
          end if
        end if
      end if
      
    end do

    if ((nvi%d(i1) .ne. 0) .and. (nvi%d(j1) .ne. 0)) then
      
      k1 = i1 + j1
        
      if (k1 > maxmass) then
        
        Allocate(nvc%d(maxmass))
        nvc%d(1:maxmass) = nvi%d(1:maxmass)
        call nvi%deinit()
        call nvi%init(maxmass*4)
        nvi%d(1:maxmass) = nvc%d(1:maxmass)
        Deallocate(nvc%d)
          
        treestart = treestart + maxmass/4
        maxmass = maxmass*4
          
        call u%deinit()
        call v%deinit()
        call utree%deinit()
        call vtree%deinit()
        call u%init(maxmass)
        call v%init(maxmass)
        call ui%init(maxmass)
        call vi%init(maxmass)
        do j = 1, maxmass
          u%d(j) = u_vec(j)
          v%d(j) = v_vec(j)
        end do
        do j = 1, maxmass
          ui%d(j) = nvi%d(j)*u%d(j)
          vi%d(j) = nvi%d(j)*v%d(j)
        end do
        utree = treebuildv(ui, treestart)
        vtree = treebuildv(vi, treestart)
      end if
        
      np = np-1
      nvi%d(i1) = nvi%d(i1)-1
      nvi%d(j1) = nvi%d(j1)-1
      nvi%d(k1) = nvi%d(k1)+1
        
      !Multiplication is faster and better than summation
      ui%d(i1) = nvi%d(i1)*u%d(i1)
      ui%d(j1) = nvi%d(j1)*u%d(j1)
      ui%d(k1) = nvi%d(k1)*u%d(k1)
      vi%d(i1) = nvi%d(i1)*v%d(i1)
      vi%d(j1) = nvi%d(j1)*v%d(j1)
      vi%d(k1) = nvi%d(k1)*v%d(k1)
        
      !Order is ok
      call treeupdatev(utree, ui, k1, treestart,u%d(k1))
      call treeupdatev(utree, ui, i1, treestart,-u%d(i1))
      call treeupdatev(utree, ui, j1, treestart,-u%d(j1))
      usum = utree%d(1) + utree%d(2) + utree%d(3) + utree%d(4)
      call treeupdatev(vtree, vi, k1, treestart,v%d(k1))
      call treeupdatev(vtree, vi, i1, treestart,-v%d(i1))
      call treeupdatev(vtree, vi, j1, treestart,-v%d(j1))
      vsum = vtree%d(1) + vtree%d(2) + vtree%d(3) + vtree%d(4)
      totalrate = usum*vsum
    end if

    if (np <= maxnp/2) then
      np = np*2
      nvi%d(:) = nvi%d(:)*2
      ui = ui*2.0d0
      vi = vi*2.0d0
      utree = treebuildv(ui, treestart)
      usum = utree%d(1) + utree%d(2) + utree%d(3) + utree%d(4)
      vtree = treebuildv(vi, treestart)
      vsum = vtree%d(1) + vtree%d(2) + vtree%d(3) + vtree%d(4)
      totalrate = usum*vsum
      nv1 = nv1/2
    end if
    
    if ((verbose > 1) .and. ((mod(step,everys) == 0) .or. (curtime >= maxtime))) then
      print *, step, curtime, np, maxmass
      print *, 'Time passed', dsecnd() - dstart
    end if
  end do
!   close(2)
  nv0 = nv1*np
  
  if (verbose > 0) then
    if (verbose == 1) then
      print *, 'Time passed', dsecnd() - dstart
    end if
    print *, 'Collisions', step
    print *, 'Number density', nv0
    print *, 'Maxmass', maxmass
    print *, 'DONE'
  end if
end

!maxmass should be a power of 4
!fastest trees.
subroutine MonteSimplePart(nv0, nvi, maxtime, verbose)
  Double precision, intent(inout) :: nv0 !Total number density; IN and OUT
  Type(IntVec) :: nvi !Numbers of particles; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  Integer(4) npsize !Size of the array of particles
  Type(IntVec) mvi !Array of particle masses
  
  Type(Vector) ui, vi !Low-rank factors
  Type(Vector) utree, vtree !Segment trees
  Double precision usum, vsum, unew, vnew !Vectors, containing sums of collision rates
  
  Integer(4) np, maxnp !Current and initial total number of particles
  Integer(4) treestart !Start index of segment tree data
  
  Integer(8) step !Current step
  Integer(8) everys !How often write to file
  Double precision nv1 !Number density per particle
  Double precision curtime !Current system time
  Double precision randi !Random number
  Double precision totalrate !Total collision rate
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  
  Integer(4) j, k, i1, j1, im, jm, k1 !Indices
  
  !Temporaries
  Type(IntVec) nvc
  
  np = sum(nvi%d)
  maxnp = np
  everys = np!/5
  
  nv1 = nv0/np

  npsize = max(16, ibset(0, 2*((bit_size(np-1) - leadz(np-1) + 1)/2)))
  
  call mvi%init(npsize)
  k = 0
  do i1 = 1, nvi%n
    do j = 1, nvi%d(i1)
      k = k + 1
      mvi%d(k) = i1
    end do
  end do
  
  call nvc%init(np/2)
  
  treestart = 0
  !Calculate geometric series
  i1 = 1
  do while (i1*4 < npsize)
    treestart = treestart + i1
    i1 = i1*4
  end do

  call ui%init(npsize)
  call vi%init(npsize)
  do j = 1, np
    ui%d(j) = u_vec(mvi%d(j))
    vi%d(j) = v_vec(mvi%d(j))
  end do
  utree = treebuildv(ui, treestart)
  usum = utree%d(1) + utree%d(2) + utree%d(3) + utree%d(4)
  vtree = treebuildv(vi, treestart)
  vsum = vtree%d(1) + vtree%d(2) + vtree%d(3) + vtree%d(4)
  totalrate = usum*vsum
  
!   open(2, file='TimeLowrankPart')
  curtime = 0.0d0
  step = 0
  dstart = dsecnd()
  do while (curtime < maxtime)
    step = step + 1

    i1 = 0
    do while (.true.)
    
      curtime = curtime + 1.0d0/totalrate/nv1
      
      randi = grnd()*usum
      i1 = treefindv(utree, ui, randi, treestart)
      randi = grnd()*vsum
      j1 = treefindv(vtree, vi, randi, treestart)
      
      if (i1 .ne. j1) then
        im = mvi%d(i1)
        jm = mvi%d(j1)
        if ((im .ne. 0) .and. (jm .ne. 0)) then
          if (grnd() <= sucrates(im, jm)) then
            exit
          end if
        end if
      end if
      
    end do
      
    k1 = im + jm
        
    np = np-1
    mvi%d(i1) = k1
    mvi%d(j1) = 0
    
    unew = u_vec(k1)
    vnew = v_vec(k1)
        
    call treeupdatev(utree, ui, i1, treestart,unew-ui%d(i1))
    call treeupdatev(utree, ui, j1, treestart,-ui%d(j1))
    usum = utree%d(1) + utree%d(2) + utree%d(3) + utree%d(4)
    call treeupdatev(vtree, vi, i1, treestart,vnew-vi%d(i1))
    call treeupdatev(vtree, vi, j1, treestart,-vi%d(j1))
    vsum = vtree%d(1) + vtree%d(2) + vtree%d(3) + vtree%d(4)
    totalrate = usum*vsum
    
    ui%d(i1) = unew
    vi%d(i1) = vnew
    ui%d(j1) = 0
    vi%d(j1) = 0

    if (np <= maxnp/2) then
      k = 1
      do j = 1, maxnp
        if (mvi%d(j) > 0) then
          nvc%d(k) = mvi%d(j)
          k = k+1
        end if
      end do
      mvi%d(1:np) = nvc%d(1:np)
      mvi%d(np+1:2*np) = nvc%d(1:np)
      do j = 1, np
        ui%d(j) = u_vec(mvi%d(j))
        vi%d(j) = v_vec(mvi%d(j))
      end do
      ui%d(np+1:2*np) = ui%d(1:np)
      vi%d(np+1:2*np) = vi%d(1:np)
      np = np*2
      utree = treebuildv(ui, treestart)
      usum = utree%d(1) + utree%d(2) + utree%d(3) + utree%d(4)
      vtree = treebuildv(vi, treestart)
      vsum = vtree%d(1) + vtree%d(2) + vtree%d(3) + vtree%d(4)
      totalrate = usum*vsum
      nv1 = nv1/2
    end if
    
    if ((verbose > 1) .and. ((mod(step,everys) == 0) .or. (curtime >= maxtime))) then
      print *, step, curtime, np
      print *, 'Time passed', dsecnd() - dstart
    end if
  end do
!   close(2)
  nv0 = nv1*np
  
  if (verbose > 0) then
    if (verbose == 1) then
      print *, 'Time passed', dsecnd() - dstart
    end if
    
    print *, 'Collisions', step
    print *, 'Number density', nv0
    print *, 'DONE'
  end if
end

subroutine MonteUpB(n0, nv, maxtime, verbose)
  Double precision, intent(in) :: n0 !Initial number density
  Integer(4), intent(in) :: nv !Initial number of monomers
  Double precision, intent(in) :: maxtime !Final time
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  
  Integer(4) nmax !Maximum particle mass
  Integer(4) np, maxnp !Current and initial total number of particles
  
  Integer(8) i !Step
  Integer(8) everys !How often write to file
  Integer(4), Allocatable :: mvi(:) !Arrays of particle masses
  Double precision nv1 !Number density per particle
  Double precision curtime !Current system time
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  Double precision lambda !Fragmentation rate
  Double precision kmax !Maximum kernel value (majorant)
  integer(4) total_d !Doublings
  
  Integer(4) j, i1, j1, i2, j2, k1 !Indices
  
  total_d = 1

  np = nv
  maxnp = np
  everys = np
  nv1 = n0/np
  nmax = 1
  lambda = 0.0d0!0.01d0

  !Allocate(mvi(32*np)) !In case fragmentation dominates
  Allocate(mvi(np))
  mvi = 1
  
  kmax = Cmax(1,1,nmax,nmax)
  
  curtime = 0.0d0
  i = 0
  dstart = dsecnd()
  do while (curtime < maxtime)
    i = i + 1
    
    do while(.true.)
    
      
      
      curtime = curtime + 2.0d0/dble(np)/dble(np-1)/nv1/(1+lambda)/kmax
      
      j2 = np+1
      i2 = floor(np*(1-grnd()))+1
      do while ((j2 > np) .or. (j2 == i2))
        j2 = floor(np*(1-grnd()))+1
      end do
      i1 = mvi(i2)
      j1 = mvi(j2)
      
      if (Cij(i1,j1) >= grnd()*kmax) then
        exit
      end if
    
    end do

    if (grnd() < lambda) then
      mvi(i2) = 1
      mvi(j2) = 1
      np = np + i1 + j1 - 2
    else
      k1 = i1 + j1
      if (k1 > nmax) then
        nmax = k1
        kmax = Cmax(1,1,nmax,nmax)
      end if
        
      mvi(min(i2,j2)) = k1
      mvi(max(i2,j2)) = mvi(np)
      mvi(np) = 1
      np = np-1
    end if

    if (np <= maxnp/2) then
      if (verbose > 1) then
        print *, 'doubling'
      end if
      !i = i * 2
      total_d = total_d*2
      do j = 1, np
        mvi(np+j) = mvi(j)
      end do
      np = np*2
      nv1 = nv1/2
    end if
  end do
  if (verbose > 0) then
    print *, 'Time passed', dsecnd() - dstart
    print *, 'Collisions', i
    print *, 'Number density', np*nv1
    print *, 'Maxmass', nmax
    print *, 'DONE'
  end if
end

subroutine MonteUpBtemp(n0, nv, temp, maxtime, verbose)
  Double precision, intent(in) :: n0 !Initial number density
  Integer(4), intent(in) :: nv !Initial number of monomers
  Double precision, intent(in) :: temp !Initial temperature
  Double precision, intent(in) :: maxtime !Final time
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  
  Integer(4) nmax, maxmass !Maximum particle mass
  Integer(4) np, maxnp !Current and initial total number of particles
  Integer(8) i !Step
  Integer(8) everys !How often write to file
  Integer(4), Allocatable :: mvi(:) !Arrays of particle masses
  Double precision nv1 !Number density per particle
  Double precision curtime !Current system time
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  Double precision lambda, apar !Kernel parameters
  Double precision kmax !Maximum kernel value (majorant)
  Double precision randi !Random number
  Integer(4) total_d !Doublings
  
  Type(Vector) temps !Temperatures
  Double precision tempmax !Maximum temperature
  Type(IntVec) nvi !Number densities
  Integer(4) j, i1, j1, i2, j2, k1 !Indices
  
  !Temporaries
  Type(IntVec) tmpivec
  Type(Vector) tmpvec
  Double precision tmp1, tmp2, tmp3
  
  total_d = 1

  np = nv
  maxnp = np
  everys = np
  nv1 = n0/np
  nmax = 1
  
  lambda = 0.4d0
  apar = 0.1d0
  maxmass = 64
  
  call temps%init(maxmass)
  call nvi%init(maxmass)
  nvi%d(1) = np
  temps%d(1) = temp
  tempmax = temp

  !Allocate(mvi(32*np))
  Allocate(mvi(np))
  mvi = 1
  
  curtime = 0.0d0
  i = 0
  dstart = dsecnd()
  do while (curtime < maxtime)
    i = i + 1
    
    do while(.true.)
    
      kmax = sqrt(pi/2)*(1 + nmax**(1.0d0/3.0d0))**2*sqrt(1+1.0d0/nmax)*sqrt(tempmax)
      curtime = curtime + 2.0d0/dble(np)/dble(np-1)/nv1/kmax
      
      j2 = np+1
      i2 = floor(np*(1-grnd()))+1
      do while ((j2 > np) .or. (j2 == i2))
        j2 = floor(np*(1-grnd()))+1
      end do
      i1 = mvi(i2)
      j1 = mvi(j2)
      
      randi = grnd()*kmax/sqrt(pi/2)
      if ((i1**(1.0d0/3.0d0)+j1**(1.0d0/3.0d0))**2*sqrt(temps%d(i1)/i1+temps%d(j1)/j1) >= randi) then
        exit
      end if
    
    end do

    k1 = i1 + j1
    if (k1 > maxmass) then
      call tmpvec%copy(temps)
      call tmpivec%copy(nvi)
      call temps%deinit()
      call nvi%deinit()
      call temps%init(maxmass*2)
      call nvi%init(maxmass*2)
      temps%d(:maxmass) = tmpvec%d(:maxmass)
      nvi%d(:maxmass) = tmpivec%d(:maxmass)
      maxmass = maxmass*2
      call tmpvec%deinit()
      call tmpivec%deinit()
    end if
    if (grnd() < BounceProb(lambda, apar, i1, j1, temps%d(i1), temps%d(j1))) then
      tmp1 = BounceMinus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), nvi%d(i1))
      tmp2 = BounceMinus(lambda, apar, j1, i1, temps%d(j1), temps%d(i1), nvi%d(j1))
      if (i1 == j1) then
        tmp2 = 2*tmp1 - temps%d(j1)
        tmp1 = tmp2
      end if
    else
      if (k1 > nmax) then
        nmax = k1
      end if
      
      tmp1 = AggloMinus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), nvi%d(i1))
      tmp2 = AggloMinus(lambda, apar, j1, i1, temps%d(j1), temps%d(i1), nvi%d(j1))
      if (i1 == j1) then
        tmp2 = max(0.0d0, 2*tmp1 - temps%d(j1))
        tmp1 = tmp2
      end if
      tmp3 = AggloPlus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), temps%d(k1), nvi%d(k1))
        
      mvi(min(i2,j2)) = k1
      mvi(max(i2,j2)) = mvi(np)
      mvi(np) = 1
      temps%d(k1) = tmp3
      nvi%d(i1) = nvi%d(i1)-1
      nvi%d(j1) = nvi%d(j1)-1
      nvi%d(k1) = nvi%d(k1)+1
      np = np-1
    end if
    temps%d(i1) = tmp1
    temps%d(j1) = tmp2
    tempmax = max(tempmax, temps%d(i1), temps%d(j1), temps%d(k1))

    if (np <= maxnp/2) then
      if (verbose > 1) then
        print *, 'doubling'
      end if
      !i = i * 2
      total_d = total_d*2
      do j = 1, np
        mvi(np+j) = mvi(j)
      end do
      np = np*2
      nvi%d(:) = nvi%d(:)*2
      tempmax = maxval(temps%d(:))
      nv1 = nv1/2
    end if
  end do
  if (verbose > 0) then
    print *, 'Time passed', dsecnd() - dstart
    print *, 'Collisions', i
    print *, 'Number density', np*nv1
    print *, 'Maxmass', nmax
    print *, 'DONE'
  end if
end

subroutine MonteInverse(n0, nvin, maxtime, verbose)
  Double precision, intent(in) :: n0 !Initial number density
  Integer(4), intent(in) :: nvin !Initial number of monomers
  Double precision, intent(in) :: maxtime !Final time
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  
  Type(Vector) nv !Number densities
  
  Integer(4) n !Maximum mass
  Integer(4) np, maxnp !Current and initial total number of particles
  Integer(8) i !Step
  Integer(8) everys !How often write to file
  Integer(4), Allocatable :: nvi(:) !Number densities
  
  Double precision nv1 !Number density per particle
  Double precision curtime !Current system time
  Double precision randi !Random number
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  
  Type(Mtrx) u, ui !Collision rates
  Type(Vector) totalvec !Sums of collision rates over rows
  Double precision totalrate !Total collision rate
  
  Integer(4) j, l, i1, j1, k1 !Indices
  
  !Temporaries
  Integer(4), Allocatable :: nvc(:)
  Double precision tau, tmp1, tmp3

  n = 4
  np = nvin
  maxnp = np
  everys = np!/10
  nv1 = n0/np
  
  call nv%init(n)
  nv%d(1) = nv1
  Allocate(nvi(n))
  nvi(:) = 0
  nvi(1) = np
  !call teta%init(n)
  !tavg = teta%d(1)

  call u%init(n, n)
  call ApproximateNo(u, n)
  call ui%init(n, n)
  call totalvec%init(n)
  do i = 1, n
    do j = 1, n
      ui%d(i,j) = nvi(j)*(nvi(i)*u%d(i,j))
      totalvec%d(i) = totalvec%d(i) + ui%d(i,j)
    end do
  end do
  totalrate = sum(totalvec%d)
  
  curtime = 0.0d0
  i = 0
  dstart = dsecnd()
  do while (curtime < maxtime)
    i = i + 1
    if ((mod(i,everys) == 0) .and. (verbose > 1)) then
      print *, i/everys, curtime, n, np
      print *, nv%d(1)
      print *, totalrate*nv1, sum(totalvec%d)*nv1, 2*np*n0
      print *, 'Time passed', dsecnd() - dstart
    end if

    j1 = 0
    do while (j1 == 0)
      !call random_number(randi)
      randi = 0
      do while(randi == 0)
        randi = grnd()
      end do
      tau = -log(randi)/totalrate/nv1*2!*2
      curtime = curtime + tau

        !call random_number(randi)
        randi = grnd()*totalrate
        
        !Line search
        i1 = n
        do j = 1, n
          randi = randi - totalvec%d(j)
          if (randi <= 0) then
            i1 = j
            exit
          end if
        end do
        randi = grnd()*totalvec%d(i1)
        do j = 1, n
          randi = randi - ui%d(i1,j)
          if (randi <= 0) then
            j1 = j
            exit
          end if
        end do
        
      if (i1 == j1) then
        !call random_number(randi)
        if (grnd() * nvi(i1) <= 1.0d0) then
          j1 = 0
        end if
      end if
    end do

    if ((nvi(i1) .ne. 0) .and. (nvi(j1) .ne. 0)) then

      !Breakage
      if (grnd() < 0.0d0) then
      
        np = np-2+i1+j1
        nvi(i1) = nvi(i1)-1
        nv%d(i1) = nvi(i1)*nv1
        nvi(j1) = nvi(j1)-1
        nv%d(j1) = nvi(j1)*nv1
        nvi(1) = nvi(1)+i1+j1
        nv%d(1) = nvi(1)*nv1
      
        !size 1
        if ((i1 .ne. 1) .and. (j1 .ne. 1)) then
          do j = 1, n
            tmp1 = 2*ui%d(1,j)
            ui%d(1,j) = (u%d(1,j)*nvi(1))*nvi(j)!*sqrt(teta%d(1)+teta%d(1))
            ui%d(j,1) = ui%d(1,j)
            totalrate = totalrate - tmp1
            totalrate = totalrate + 2*ui%d(1,j)
            totalvec%d(j) = totalvec%d(j) + ui%d(1,j) - tmp1/2
            totalvec%d(1) = totalvec%d(1) + ui%d(1,j) - tmp1/2
            if (j == 1) then
              totalrate = totalrate + tmp1/2 - ui%d(1,j)
              totalvec%d(1) = totalvec%d(1) - ui%d(1,j) + tmp1/2
            end if
          end do
        end if
        
        k1 = 0
        
      else
      
        k1 = i1 + j1
        if (k1 > 256000000) then
          exit
        end if
        if (k1 > n) then
          Allocate(nvc(n))
          nvc(1:n) = nvi(1:n)
          Deallocate(nvi)
          Allocate(nvi(n*2))
          nvi(1:n) = nvc(1:n)
          nvi(n+1:n*2) = 0
          Deallocate(nvc)
          call nv%deinit()
          call nv%init(n*2)
          nv%d(1:n) = nvi(1:n)*nv1
          call ApproximateNo(u, n*2)
          call ui%deinit()
          call ui%init(n*2,n*2)
          call totalvec%deinit()
          call totalvec%init(n*2)
          do j = 1, n
            do l = 1, n
              ui%d(l,j) = nvi(l)*(nvi(j)*u%d(l,j))!*sqrt(teta%d(l)+teta%d(j))
              totalvec%d(j) = totalvec%d(j) + ui%d(l,j)
            end do
          end do
          if (verbose > 1) then
            print *, 'N updated. Totalrate (new, old):', totalrate, sum(totalvec%d)
            print *, 'Difference:', totalrate-sum(totalvec%d)
          end if
          totalrate = sum(totalvec%d)
          n = n*2
        end if
        
        np = np-1
        nvi(i1) = nvi(i1)-1
        nv%d(i1) = nvi(i1)*nv1
        nvi(j1) = nvi(j1)-1
        nv%d(j1) = nvi(j1)*nv1
        nvi(k1) = nvi(k1)+1
        nv%d(k1) = nvi(k1)*nv1
        do j = 1, n
          tmp3 = 2*ui%d(k1,j)
          ui%d(k1,j) = (u%d(k1,j)*nvi(k1))*nvi(j)
          ui%d(j,k1) = ui%d(k1,j)
          totalrate = totalrate - tmp3
          totalrate = totalrate + 2*ui%d(k1,j)
          totalvec%d(j) = totalvec%d(j) + ui%d(k1,j) - tmp3/2
          totalvec%d(k1) = totalvec%d(k1) + ui%d(k1,j) - tmp3/2
          if (j == k1) then
            totalrate = totalrate + tmp3/2 - ui%d(k1,j)
            totalvec%d(k1) = totalvec%d(k1) - ui%d(k1,j) + tmp3/2
          end if
        end do
      end if
      do j = 1, n
        tmp1 = 2*ui%d(i1,j)
        ui%d(i1,j) = (u%d(i1,j)*nvi(i1))*nvi(j)
        ui%d(j,i1) = ui%d(i1,j)
        totalrate = totalrate - tmp1
        totalrate = totalrate + 2*ui%d(i1,j)
        totalvec%d(j) = totalvec%d(j) + ui%d(i1,j) - tmp1/2
        totalvec%d(i1) = totalvec%d(i1) + ui%d(i1,j) - tmp1/2
        if (j == i1) then
          totalrate = totalrate + tmp1/2 - ui%d(i1,j)
          totalvec%d(i1) = totalvec%d(i1) - ui%d(i1,j) + tmp1/2
        end if
      end do
      if (i1 .ne. j1) then
        do j = 1, n
          tmp1 = 2*ui%d(j1,j)
          ui%d(j1,j) = (u%d(j1,j)*nvi(j1))*nvi(j)
          ui%d(j,j1) = ui%d(j1,j)
          totalrate = totalrate - tmp1
          totalrate = totalrate + 2*ui%d(j1,j)
          totalvec%d(j) = totalvec%d(j) + ui%d(j1,j) - tmp1/2
          totalvec%d(j1) = totalvec%d(j1) + ui%d(j1,j) - tmp1/2
          if (j == j1) then
            totalrate = totalrate + tmp1/2 - ui%d(j1,j)
            totalvec%d(j1) = totalvec%d(j1) - ui%d(j1,j) + tmp1/2
          end if
        end do
      end if
    end if

    if (np <= maxnp/2) then
      np = np*2
      nvi = nvi*2
      ui = ui*4.0d0
      totalrate = totalrate*4.0d0
      totalvec = totalvec*4.0d0
      nv1 = nv1/2
    end if
  end do
  if (verbose > 0) then
    print *, 'Time passed', dsecnd() - dstart
    print *, 'Collisions', i
    print *, 'Number density', np*nv1
    print *, 'Maxmass', n
    print *, 'DONE'
  end if
end

!Returns kernel without any low-rank approximation
subroutine ApproximateNo(u, n)
  Type(Mtrx) :: u
  Integer(4) :: n
  Integer(4) i, j

  call u%deinit()
  call u%init(n, n)
  do i = 1, n
    do j = 1, n
      u%d(i,j) = Cij(i,j)
    end do
  end do
end

subroutine MonteWagnerRanked(r, nv0, nvi, temps, maxtime, lambda, apar, verbose)
  Integer(4), intent(in) :: r !Majorant rank
  Double precision, intent(inout) :: nv0 !Total number density; IN and OUT
  Type(IntVec) :: nvi !Numbers of particles; IN and OUT
  Type(Vector) :: temps !Temperatures of each size; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Double precision, intent(in) :: lambda, apar !Kernel parameters
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  
  Integer(4) maxmass !Maximum cluster mass in the system
  Integer(4) np, maxnp !Current and initial total number of particles
  Integer(8) step !Current step
  Integer(8) everys !How often write to file
  Integer(8) everys2 !How often apply thermostat
  Integer(4) curr !Iteration over ranks
  Double precision nv1 !Number density per particle
  Double precision tavg !Current average temperature
  Double precision timeold, curtime !Current and previous thermostat system time
  Double precision randi !Random number
  Double precision totalrate !Total collision rate
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  Type(Vector) totali, totalj, utmp, vtmp, totalvec !Vectors, containing sums of collision rates
  Double precision gsmall, Tmolgas !Molecular gas parameters
  
  Integer(4), parameter :: maxmax = 32 !coef**maxmax is the maximum possible mass.
  Type(IntVec) partof(maxmax) !Particle bins
  Type(Mtrx) ratesi, ratesj, u, v !Collision rate factors
  Integer(4) maxpartof !Current maximum possible index i in partof(i)
  Integer(8) maxpartofmass !Current maximum possible mass
  Integer(4) coef !Mass multiplier between bins
  Type(IntVec) maxmasses !Maximum mass in each bin
  
  Integer(4) i3min, j, i1, j1, k1, i2, j2, i3, j3, k3 !Indices
  
  !Temporaries
  Double precision tau, tmp1, tmp2, tmp3
  Type(Vector) tmpvec
  Type(IntVec) nvc
  
  coef = 2

  maxmass = nvi%n
  
  !Molecular gas effect
  gsmall = -1.0d0/3.0d0 !(-4/3)
  Tmolgas = 0.0d0!1.0d0
  
  np = sum(nvi%d)
  maxnp = np
  everys = np!/5
  everys2 = np!/100
  
  nv1 = nv0/np
  
  maxpartof = floor(log(maxmass-0.5d0)/log(dble(coef))) + 1
  maxpartofmass = coef**maxpartof-1
  do j = 1, maxpartof
    call partof(j)%init(np)
  end do
  call maxmasses%init(maxmax)
  call ratesi%init(r,maxmax)
  call ratesj%init(r,maxmax)
  call u%init(r,maxmax)
  call v%init(r,maxmax)
  call totali%init(r)
  call totalj%init(r)
  call utmp%init(r)
  call vtmp%init(r)
  call totalvec%init(r)
  
  i3min = 1
  
  j1 = 1
  do j = 1, maxmass
    if (coef**j1 <= j) then
      j1 = j1 + 1
    end if
    do i1 = 1, nvi%d(j)
      maxmasses%d(j1) = maxmasses%d(j1) + 1
      partof(j1)%d(maxmasses%d(j1)) = j
    end do
    call Cuv(j,temps%d(j),utmp,vtmp)
    
    u%d(1:r,j1) = max(u%d(1:r,j1),utmp%d(1:r))
    v%d(1:r,j1) = max(v%d(1:r,j1),vtmp%d(1:r))
  end do
  
  do j = 1, maxpartof
    ratesi%d(1:r,j) = maxmasses%d(j)*u%d(1:r,j)
  end do
  do j = 1, maxpartof
    ratesj%d(1:r,j) = maxmasses%d(j)*v%d(1:r,j)
  end do
  
  totali%d = sum(ratesi%d,2)
  totalj%d = sum(ratesj%d,2)
  totalvec%d(1:r) = totali%d(1:r)*totalj%d(1:r)
  totalrate = sum(totalvec%d)
  
  timeold = curtime
  
!   open(2, file='TimeWagnerTemp')
  curtime = 0.0d0
  step = 0
  dstart = dsecnd()
  do while (curtime < maxtime)
    step = step + 1

    i1 = 0
    do while (1 > 0)
    
      tau = 1.0d0/totalrate/nv1
      curtime = curtime + tau
      
      randi = grnd()*totalrate
      do curr = 1, r-1
        randi = randi - totalvec%d(curr)
        if (randi <= 0) then
          exit
        end if
      end do
      
        i3 = maxpartof+1
        do while (i3 > maxpartof)
          randi = grnd()*totali%d(curr)
          do i3 = i3min, maxpartof
            randi = randi - ratesi%d(curr,i3)
            if (randi <= 0) then
              exit
            end if
          end do
        end do
        
        j3 = maxpartof+1
        do while (j3 > maxpartof)
          randi = grnd()*totalj%d(curr)
          do j3 = i3min, maxpartof
            randi = randi - ratesj%d(curr,j3)
            if (randi <= 0) then
              exit
            end if
          end do
        end do
        
        i2 = maxmasses%d(i3)+1
        j2 = maxmasses%d(j3)+1
        do while (i2 > maxmasses%d(i3))
          i2 = floor(grnd()*maxmasses%d(i3))+1
        end do
        do while (j2 > maxmasses%d(j3))
          j2 = floor(grnd()*maxmasses%d(j3))+1
        end do
        
        i1 = partof(i3)%d(i2)
        j1 = partof(j3)%d(j2)
      
      if (grnd() <= CijTemp(i1,j1,temps%d(i1),temps%d(j1))/(sum(u%d(1:r,i3)*v%d(1:r,j3))+sum(u%d(1:r,j3)*v%d(1:r,i3)))) then
        if ((i2 .ne. j2) .or. (i3 .ne. j3)) then
          exit
        end if
      end if
      
    end do

    if ((temps%d(i1) + temps%d(j1) > 0) .and. (nvi%d(i1) .ne. 0) .and. (nvi%d(j1) .ne. 0)) then
      if (grnd() < BounceProb(lambda, apar, i1, j1, temps%d(i1), temps%d(j1))) then
      
        tmp1 = BounceMinus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), nvi%d(i1))
        tmp2 = BounceMinus(lambda, apar, j1, i1, temps%d(j1), temps%d(i1), nvi%d(j1))
        if (i1 == j1) then
          tmp2 = 2*tmp1 - temps%d(j1)
          tmp1 = tmp2
        end if
        
      else
      
        k1 = i1 + j1
        
        if (k1 > maxmass) then
        
          call tmpvec%copy(temps)
          call temps%deinit()
          call temps%init(maxmass*2)
          temps%d(1:maxmass) = tmpvec%d(1:maxmass)
          temps%d(maxmass+1:maxmass*2) = 1
          call tmpvec%deinit()
          Allocate(nvc%d(maxmass))
          nvc%d(1:maxmass) = nvi%d(1:maxmass)
          call nvi%deinit()
          call nvi%init(maxmass*2)
          nvi%d(1:maxmass) = nvc%d(1:maxmass)
          Deallocate(nvc%d)
          
          maxmass = maxmass*2
          
        end if
        
        if (k1 > maxpartofmass) then
          maxpartof = maxpartof + 1
          maxpartofmass = coef**maxpartof - 1
          call partof(maxpartof)%init(maxnp)
          do while (maxmasses%d(i3min) == 0)
            i3min = i3min+1
          end do
        end if
        
        tmp1 = AggloMinus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), nvi%d(i1))
        tmp2 = AggloMinus(lambda, apar, j1, i1, temps%d(j1), temps%d(i1), nvi%d(j1))
        if (i1 == j1) then
          tmp2 = 2*tmp1 - temps%d(j1)
          tmp1 = tmp2
        end if
        tmp3 = AggloPlus(lambda, apar, i1, j1, temps%d(i1), temps%d(j1), temps%d(k1), nvi%d(k1))
        
        np = np-1
        nvi%d(i1) = nvi%d(i1)-1
        nvi%d(j1) = nvi%d(j1)-1
        nvi%d(k1) = nvi%d(k1)+1
        
        temps%d(k1) = tmp3
        
        k3 = max(i3,j3)
        
        if (i2 < j2) then
          partof(j3)%d(j2) = partof(j3)%d(maxmasses%d(j3))
          maxmasses%d(j3) = maxmasses%d(j3)-1
          partof(i3)%d(i2) = partof(i3)%d(maxmasses%d(i3))
          maxmasses%d(i3) = maxmasses%d(i3)-1
        else
          partof(i3)%d(i2) = partof(i3)%d(maxmasses%d(i3))
          maxmasses%d(i3) = maxmasses%d(i3)-1
          partof(j3)%d(j2) = partof(j3)%d(maxmasses%d(j3))
          maxmasses%d(j3) = maxmasses%d(j3)-1
        end if
        if (k1 > coef**k3 - 1) then
          k3 = k3 + 1
        end if
        maxmasses%d(k3) = maxmasses%d(k3)+1
        partof(k3)%d(maxmasses%d(k3)) = k1
        
        call Cuv(k1,temps%d(k1),utmp,vtmp)
        u%d(1:r,k3) = max(u%d(1:r,k3),utmp%d(1:r))
        v%d(1:r,k3) = max(v%d(1:r,k3),vtmp%d(1:r))
        
        totali%d(1:r) = totali%d(1:r) - ratesi%d(1:r,k3)
        ratesi%d(1:r,k3) = maxmasses%d(k3)*u%d(1:r,k3)
        totali%d(1:r) = totali%d(1:r) + ratesi%d(1:r,k3)
  
        totalj%d(1:r) = totalj%d(1:r) - ratesj%d(1:r,k3)
        ratesj%d(1:r,k3) = maxmasses%d(k3)*v%d(1:r,k3)
        totalj%d(1:r) = totalj%d(1:r) + ratesj%d(1:r,k3)
        
      end if
      
      temps%d(i1) = tmp1
      temps%d(j1) = tmp2
      
      call Cuv(i1,temps%d(i1),utmp,vtmp)
      u%d(1:r,i3) = max(u%d(1:r,i3),utmp%d(1:r))
      v%d(1:r,i3) = max(v%d(1:r,i3),vtmp%d(1:r))
      call Cuv(j1,temps%d(j1),utmp,vtmp)
      u%d(1:r,j3) = max(u%d(1:r,j3),utmp%d(1:r))
      v%d(1:r,j3) = max(v%d(1:r,j3),vtmp%d(1:r))
      
      totali%d(1:r) = totali%d(1:r) - ratesi%d(1:r,i3)
      ratesi%d(1:r,i3) = maxmasses%d(i3)*u%d(1:r,i3)
      totali%d(1:r) = totali%d(1:r) + ratesi%d(1:r,i3)
      totali%d(1:r) = totali%d(1:r) - ratesi%d(1:r,j3)
      ratesi%d(1:r,j3) = maxmasses%d(j3)*u%d(1:r,j3)
      totali%d(1:r) = totali%d(1:r) + ratesi%d(1:r,j3)
  
      totalj%d(1:r) = totalj%d(1:r) - ratesj%d(1:r,i3)
      ratesj%d(1:r,i3) = maxmasses%d(i3)*v%d(1:r,i3)
      totalj%d(1:r) = totalj%d(1:r) + ratesj%d(1:r,i3)
      totalj%d(1:r) = totalj%d(1:r) - ratesj%d(1:r,j3)
      ratesj%d(1:r,j3) = maxmasses%d(j3)*v%d(1:r,j3)
      totalj%d(1:r) = totalj%d(1:r) + ratesj%d(1:r,j3)
      
      totalvec%d(1:r) = totali%d(1:r)*totalj%d(1:r)
      totalrate = sum(totalvec%d)
      
    end if

    if (np <= maxnp/2) then
      np = np*2
      nvi%d(:) = nvi%d(:)*2
      do j = 1, maxpartof
        do j1 = 1, maxmasses%d(j)
          partof(j)%d(j1+maxmasses%d(j)) = partof(j)%d(j1)
        end do
        maxmasses%d(j) = 2*maxmasses%d(j)
      end do
      
      u%d(1:r,1:maxpartof) = 0
      v%d(1:r,1:maxpartof) = 0
      j1 = 1
      do j = 1, maxmass
        if (coef**j1 <= j) then
          j1 = j1 + 1
        end if
        if (nvi%d(j) > 0) then
          call Cuv(j,temps%d(j),utmp,vtmp)
          u%d(1:r,j1) = max(u%d(1:r,j1),utmp%d(1:r))
          v%d(1:r,j1) = max(v%d(1:r,j1),vtmp%d(1:r))
        end if
      end do
  
      do j = 1, maxpartof
        ratesi%d(1:r,j) = maxmasses%d(j)*u%d(1:r,j)
      end do
      do j = 1, maxpartof
        ratesj%d(1:r,j) = maxmasses%d(j)*v%d(1:r,j)
      end do
  
      totali%d = sum(ratesi%d,2)
      totalj%d = sum(ratesj%d,2)
      totalvec%d(1:r) = totali%d(1:r)*totalj%d(1:r)
      totalrate = sum(totalvec%d)

      nv1 = nv1/2
    end if
    
    if ((mod(step,everys2) == 0) .or. (curtime >= maxtime)) then
    
      tavg = 0.0d0
      do j = 1, maxmass
        tavg = tavg + temps%d(j)*nvi%d(j)
      end do
      tavg = tavg/np
      
      timeold = curtime
      
      if (verbose > 1) then
!         write(2, *) curtime, temps%d(1)!, temps%d(3), temps%d(5), temps%d(10), temps%d(20), tavg
      end if
    end if
    if ((verbose > 1) .and. ((mod(step,everys) == 0) .or. (curtime >= maxtime))) then
      print *, step, curtime, np, maxmass
      print *, 'Temperatures', temps%d(1), temps%d(2), tavg
      print *, 'Time passed', dsecnd() - dstart
    end if
  end do
!   close(2)
  nv0 = nv1*np
  
  if (verbose > 0) then
    if (verbose == 1) then
      print *, 'Time passed', dsecnd() - dstart
    end if
    print *, 'Collisions', step
  
    print *, 'Number density', nv0
    print *, 'Tavg', tavg
    print *, 'Maxmass', maxmass
    print *, 'DONE'
  end if
end

subroutine initWalker(small, large, prob, alias, probavg)
  Type(IntVec) :: small, large
  Type(Vector) :: prob
  Type(IntVec) :: alias
  Double precision :: probavg

  Integer(4) l, s, j, k
  
  l = 0
  s = 0
  do j = 1, prob%n
    if (prob%d(j) > probavg) then
      l = l + 1
      large%d(l) = j
    else
      s = s + 1
      small%d(s) = j
    end if
  end do
  do while ((s .ne. 0) .and. (l .ne. 0))
    j = small%d(s)
    k = large%d(l)
    alias%d(j) = k
    prob%d(k) = prob%d(k) + (prob%d(j) - probavg)
    if (prob%d(k) >= probavg) then
      s = s - 1
    else
      l = l - 1
      small%d(s) = k
    end if
  end do
  prob%d(small%d(:s)) = probavg
end subroutine

function getWalker(prob, alias, probavg) Result(res)
  Type(Vector), intent(in) :: prob
  Type(IntVec), intent(in) :: alias
  Double precision, intent(in) :: probavg
  Integer(4) :: res
  
  Integer(4) i
  
  i = floor(prob%n*(1-grnd()))+1
  if (grnd()*probavg <= prob%d(i)) then
    res = i
  else
    res = alias%d(i)
  end if
end function

subroutine MonteWalker(nv0, nvi, maxtime, verbose)
  Double precision, intent(inout) :: nv0 !Total number density; IN and OUT
  Type(IntVec) :: nvi !Numbers of particles; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  
  Integer(4) np, maxnp !Current and initial total number of particles
  Integer(4) maxmass !Maximum cluster mass in the system
  Integer(8) step !Current step
  Integer(8) everys !How often write to file
  Double precision nv1 !Number density per particle
  Double precision curtime !Current system time
  Double precision totalrate !Total collision rate
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  
  Integer(4), parameter :: maxmax = 128 !1/coef**maxmax is the maximum possible mass.
  Type(IntVec), allocatable :: partof(:) !Particle bins
  Integer(4) maxpartof !Current maximum possible index i in partof(i)
  Integer(4) maxold !Old maxpartof, needed for update of maxpartof
  Double precision coef !Mass multiplier between bins
  Double precision cupdate !When cupdate ratio is reached, update is run
  Type(IntVec) numbers, maxnumbers, minnumbers !Number of particles in bins, upper and lower bounds
  Double precision pcoef !Probability coefficient. Depends on whether the same bin is selected or not.
  Type(IntVec) :: lefts, rights !Bin mass bounds
  Logical toupdate !Set to True when update is needed
  Type(IntVec) masspartof !Pointers from masses to bins
  
  Integer(4) totupdates !Total probability distribution updates
  
  Integer(4) i, j, k, i1, j1, i2, j2, k2, i3, j3 !Indices
  
  !Walker method variables
  Type(IntVec) :: small, large
  Type(Vector) :: prob0, prob
  Type(IntVec) :: alias
  Double precision probavg
  
  !Temporaries
  Type(IntVec) tmpvec
  
  totupdates = 0
  
  maxmass = nvi%n
  if (verbose > 0) then
    print *, 'MaxMass', maxmass
  end if
  
  np = sum(nvi%d)
  maxnp = np
  everys = np!/5
  
  coef = min(0.99d0,max(0.2d0,4*np**(-1.0d0/6.0d0))) !0.25 for np = 10^7, 0.5 for np = 10^5
  cupdate = 1.0d0 + np**(-1.0d0/6.0d0)/3 !tested on np = 10^7
  
  Allocate(partof(maxmax))
  
  nv1 = nv0/np
  call lefts%init(maxmax) !decrease and restart them every time?
  call rights%init(maxmax)
  call numbers%init(maxmax)
  call minnumbers%init(maxmax)
  call maxnumbers%init(maxmax)
  j = 1
  lefts%d(1) = 1
  rights%d(1) = floor(1.0d0+coef)
  do while (rights%d(j) < maxmass)
    j = j+1
    lefts%d(j) = rights%d(j-1)+1
    rights%d(j) = floor(lefts%d(j)*(1.0d0 + coef))
  end do
  maxpartof = j
  
  call tmpvec%init(np)
  do j = 1, maxpartof
    do j1 = lefts%d(j), min(rights%d(j), maxmass)
      do i1 = numbers%d(j)+1, numbers%d(j) + nvi%d(j)
        tmpvec%d(i1) = j1
      end do
      numbers%d(j) = numbers%d(j) + nvi%d(j1)
    end do
    minnumbers%d(j) = floor(numbers%d(j)/cupdate)+1
    maxnumbers%d(j) = floor(cupdate*numbers%d(j))
    call partof(j)%init(min(np,4*(maxnumbers%d(j)+1)))
    partof(j)%d(1:numbers%d(j)) = tmpvec%d(1:numbers%d(j))
  end do
  maxmass = rights%d(maxpartof)
  call masspartof%init(2*maxmass)
  do j = 1, maxpartof
    do j1 = lefts%d(j), rights%d(j)
      masspartof%d(j1) = j
    end do
  end do
  
  call small%init(maxpartof**2)
  call large%init(maxpartof**2)
  call prob0%init(maxpartof**2)
  call prob%init(maxpartof**2)
  call alias%init(maxpartof**2)
  
  do i = 1, maxpartof
    do j = 1, maxpartof
      prob0%d((i-1)*maxpartof + j) = Cmax(lefts%d(i), lefts%d(j), rights%d(i), rights%d(j))*maxnumbers%d(i)*maxnumbers%d(j)
    end do
  end do
  !Try without rescale, but full reject k == l?
  do i = 1, maxpartof
    prob0%d((i-1)*maxpartof + i) = Cmax(lefts%d(i), lefts%d(i), rights%d(i), rights%d(i)) &
                                   *maxnumbers%d(i)*(maxnumbers%d(i)-1)
  end do
  
  totalrate = sum(prob0%d)
  probavg = totalrate/maxpartof**2
  call prob%copy(prob0) !or prob=prob0(:)?
  call initWalker(small, large, prob, alias, probavg)
  
!   open(2, file='TimeWalker')
  curtime = 0.0d0
  step = 0
  dstart = dsecnd()
  do while (.true.)
    step = step + 1
    
    do while (.true.)
    
      curtime = curtime + 2.0d0/totalrate/nv1
      
      i3 = getWalker(prob, alias, probavg)
      j2 = (i3-1)/maxpartof
      i2 = i3 - j2*maxpartof
      j2 = j2 + 1
      
      i1 = floor(numbers%d(i2)*(1-grnd()))+1
      if (i2 == j2) then
        j1 = i1 + floor((numbers%d(i2)-1)*(1-grnd()))+1
        if (j1 > numbers%d(i2)) then
          j1 = j1 - numbers%d(i2)
        end if
        if (j1 > i1) then
          j3 = j1
          j1 = i1
          i1 = j3
          j3 = j2
          j2 = i2
          i2 = j3
        end if
        pcoef = dble(numbers%d(i2))*(numbers%d(i2)-1)
      else
        j1 = floor(numbers%d(j2)*(1-grnd()))+1
        pcoef = dble(numbers%d(i2))*numbers%d(j2)
      end if
      
      i = partof(i2)%d(i1)
      j = partof(j2)%d(j1)
      
      if (grnd()*prob0%d(i3) <= Cij(i,j)*pcoef ) then
        exit
      end if      

    end do
    
    if (((mod(step,everys) == 0) .or. (curtime >= maxtime)) .and. (verbose == 2)) then
      print *, step, curtime, np, maxmass
      print *, 'Time passed', dsecnd() - dstart
    end if
    if (curtime >= maxtime) then
      exit
    end if

    !if ((temps%d(i1) + temps%d(j1) > 0) .and. (nvi%d(i1) .ne. 0) .and. (nvi%d(j1) .ne. 0)) then
      
        k = i + j
        
        if (k > maxmass) then
        
          maxmass = k
          
          maxold = maxpartof
          do while (rights%d(maxpartof) < maxmass)
            maxpartof = maxpartof+1
            lefts%d(maxpartof) = rights%d(maxpartof-1)+1
            rights%d(maxpartof) = floor(lefts%d(maxpartof)*(1.0d0 + coef))
          end do
  
          do j3 = maxold+1, maxpartof
            call partof(j3)%init(4)
          end do
          maxmass = rights%d(maxpartof)
          
          if (masspartof%n < maxmass) then
            call masspartof%deinit()
            call masspartof%init(2*maxmass)
            do i3 = 1, maxold
              do j3 = lefts%d(i3), rights%d(i3)
                masspartof%d(j3) = i3
              end do
            end do
          end if
          do i3 = maxold+1, maxpartof
            do j3 = lefts%d(i3), rights%d(i3)
              masspartof%d(j3) = i3
            end do
          end do
  
          call small%deinit()
          call large%deinit()
          call prob0%deinit()
          call prob%deinit()
          call alias%deinit()
          call small%init(maxpartof**2)
          call large%init(maxpartof**2)
          call prob0%init(maxpartof**2)
          call prob%init(maxpartof**2)
          call alias%init(maxpartof**2)
          
          k2 = maxpartof
        else
          k2 = masspartof%d(k)
        end if
        
        np = np-1
        
        
        partof(i2)%d(i1) = partof(i2)%d(numbers%d(i2))
        numbers%d(i2) = numbers%d(i2) - 1
        partof(j2)%d(j1) = partof(j2)%d(numbers%d(j2))
        numbers%d(j2) = numbers%d(j2) - 1
        numbers%d(k2) = numbers%d(k2) + 1
        partof(k2)%d(numbers%d(k2)) = k
        
        if (np <= maxnp/2) then
          np = np*2
          nv1 = nv1/2
          do j = 1, maxpartof
            partof(j)%d(1+numbers%d(j):2*numbers%d(j)) = partof(j)%d(1:numbers%d(j))
            numbers%d(j) = 2*numbers%d(j)
          end do
          toupdate = .true.
        end if
        
        if ((numbers%d(i2) < minnumbers%d(i2)) .or. (numbers%d(j2) < minnumbers%d(j2)) .or. (numbers%d(k2) > maxnumbers%d(k2))) then
          toupdate = .true.
        end if
        
        if (toupdate) then
          totupdates = totupdates + 1
          toupdate = .false.
          do j = 1, maxpartof
            minnumbers%d(j) = floor(numbers%d(j)/cupdate)+1
            maxnumbers%d(j) = floor(cupdate*numbers%d(j))
            if (min(maxnp,2*(maxnumbers%d(j)+1)) > partof(j)%n) then
              tmpvec%d(1:numbers%d(j)) = partof(j)%d(1:numbers%d(j))
              call partof(j)%deinit()
              call partof(j)%init(min(maxnp,4*(maxnumbers%d(j)+1)))
              partof(j)%d(1:numbers%d(j)) = tmpvec%d(1:numbers%d(j))
            end if
            if (8*(maxnumbers%d(j)+1) < partof(j)%n) then
              tmpvec%d(1:numbers%d(j)) = partof(j)%d(1:numbers%d(j))
              call partof(j)%deinit()
              call partof(j)%init(min(maxnp,4*(maxnumbers%d(j)+1)))
              partof(j)%d(1:numbers%d(j)) = tmpvec%d(1:numbers%d(j))
            end if
          end do
          
          do i = 1, maxpartof
            do j = 1, maxpartof
              prob0%d((i-1)*maxpartof + j) = Cmax(lefts%d(i), lefts%d(j), rights%d(i), rights%d(j))*maxnumbers%d(i)*maxnumbers%d(j)
            end do
          end do
          !Try without rescale, but full reject k == l?
          do i = 1, maxpartof
            prob0%d((i-1)*maxpartof + i) = Cmax(lefts%d(i), lefts%d(i), rights%d(i), rights%d(i)) &
                                           *maxnumbers%d(i)*(maxnumbers%d(i)-1)
          end do
  
          totalrate = sum(prob0%d)
          probavg = totalrate/maxpartof**2
          call prob%copy(prob0) !or prob=prob0(:)?
          call initWalker(small, large, prob, alias, probavg)
        end if
      
    !end if
    
  end do
!   close(2)
  nv0 = nv1*np
  
  if (verbose > 0) then
    print *, 'Time passed', dsecnd() - dstart
    print *, 'Collisions', step-1
    print *, 'Number density', nv0
    print *, 'Maxmass, Parts, Updates', maxmass, maxpartof**2, totupdates
    print *, 'DONE'
  end if
end

subroutine MonteWagnerSimple(nv0, nvi, maxtime, verbose)
  Double precision, intent(inout) :: nv0 !Total number density; IN and OUT
  Type(IntVec) :: nvi !Numbers of particles; IN and OUT
  Double precision, intent(in) :: maxtime !Max laboratory time
  Integer(4), intent(in) :: verbose !0 - only in files; 1 - start and finish; 2 - every everys steps
  Integer(4) maxmass !Maximum cluster mass in the system
  Integer(4) np, maxnp !Current and initial total number of particles
  Integer(8) step !Current step
  Integer(8) everys !How often write to file
  Double precision nv1 !Number density per particle
  Double precision curtime !Current system time
  Double precision randi !Random number
  Double precision totalrate !Total collision rate
  Double precision dstart !Initial time
  Double precision dsecnd !LAPACK procedure for time measurement
  Integer(4), parameter :: maxmax = 32 !coef**maxmax is the maximum possible mass.
  Type(IntVec) partof(maxmax) !Particle bins
  Type(Vector) ratesi, ratesj, u, v !Collision rate factors
  Integer(4) maxpartof !Current maximum possible index i in partof(i)
  Integer(8) maxpartofmass !Current maximum possible mass
  Integer(4) coef !Mass multiplier between bins
  Type(IntVec) maxmasses !Maximum mass in each bin
  Double precision totali, totalj !Sums of collision rates
  Integer(4) i3min, j, i1, j1, k1, i2, j2, i3, j3, k3 !Indices
  
  !Temporaries
  Double precision tau
  Type(IntVec) nvc
  
  coef = 2

  maxmass = nvi%n
  
  np = sum(nvi%d)
  maxnp = np
  everys = np!/5
  
  nv1 = nv0/np
  
  maxpartof = floor(log(maxmass-0.5d0)/log(dble(coef))) + 1
  maxpartofmass = coef**maxpartof-1
  do j = 1, maxpartof
    call partof(j)%init(np)
  end do
  call maxmasses%init(maxmax)
  call ratesi%init(maxmax)
  call ratesj%init(maxmax)
  call u%init(maxmax)
  call v%init(maxmax)
  
  i3min = 1
  
  j1 = 1
  do j = 1, maxmass
    if (coef**j1 <= j) then
      j1 = j1 + 1
    end if
    do i1 = 1, nvi%d(j)
      maxmasses%d(j1) = maxmasses%d(j1) + 1
      partof(j1)%d(maxmasses%d(j1)) = j
    end do
    u%d(j1) = max(u%d(j1),u_vec(j))
    v%d(j1) = max(v%d(j1),v_vec(j))
  end do
  
  do j = 1, maxpartof
    ratesi%d(j) = maxmasses%d(j)*u%d(j)
  end do
  do j = 1, maxpartof
    ratesj%d(j) = maxmasses%d(j)*v%d(j)
  end do
  
  totali = sum(ratesi%d)
  totalj = sum(ratesj%d)
  totalrate = totali*totalj
  
!   open(2, file='TimeWalker')
  curtime = 0.0d0
  step = 0
  dstart = dsecnd()
  do while (curtime < maxtime)
    step = step + 1

    i1 = 0
    do while (1 > 0)
    
      tau = 1.0d0/totalrate/nv1
      curtime = curtime + tau
      
        i3 = maxpartof+1
        do while (i3 > maxpartof)
          randi = grnd()*totali
          do i3 = i3min, maxpartof
            randi = randi - ratesi%d(i3)
            if (randi <= 0) then
              exit
            end if
          end do
        end do
        
        j3 = maxpartof+1
        do while (j3 > maxpartof)
          randi = grnd()*totalj
          do j3 = i3min, maxpartof
            randi = randi - ratesj%d(j3)
            if (randi <= 0) then
              exit
            end if
          end do
        end do
        
        i2 = maxmasses%d(i3)+1
        j2 = maxmasses%d(j3)+1
        do while (i2 > maxmasses%d(i3))
          i2 = floor(grnd()*maxmasses%d(i3))+1
        end do
        do while (j2 > maxmasses%d(j3))
          j2 = floor(grnd()*maxmasses%d(j3))+1
        end do
        
        i1 = partof(i3)%d(i2)
        j1 = partof(j3)%d(j2)
      
      if (grnd() <= Cij(i1,j1)/(u%d(i3)*v%d(j3)+u%d(j3)*v%d(i3))) then
        if ((i2 .ne. j2) .or. (i3 .ne. j3)) then
          exit
        end if
      end if
      
    end do

    if ((nvi%d(i1) .ne. 0) .and. (nvi%d(j1) .ne. 0)) then
      
        k1 = i1 + j1
        
        if (k1 > maxmass) then
        
          Allocate(nvc%d(maxmass))
          nvc%d(1:maxmass) = nvi%d(1:maxmass)
          call nvi%deinit()
          call nvi%init(maxmass*2)
          nvi%d(1:maxmass) = nvc%d(1:maxmass)
          Deallocate(nvc%d)
          
          maxmass = maxmass*2
          
        end if
        
        if (k1 > maxpartofmass) then
          maxpartof = maxpartof + 1
          maxpartofmass = coef**maxpartof - 1
          call partof(maxpartof)%init(maxnp)
          do while (maxmasses%d(i3min) == 0)
            i3min = i3min+1
          end do
        end if
        
        np = np-1
        nvi%d(i1) = nvi%d(i1)-1
        nvi%d(j1) = nvi%d(j1)-1
        nvi%d(k1) = nvi%d(k1)+1
        
        k3 = max(i3,j3)
        
        if (i2 < j2) then
          partof(j3)%d(j2) = partof(j3)%d(maxmasses%d(j3))
          maxmasses%d(j3) = maxmasses%d(j3)-1
          partof(i3)%d(i2) = partof(i3)%d(maxmasses%d(i3))
          maxmasses%d(i3) = maxmasses%d(i3)-1
        else
          partof(i3)%d(i2) = partof(i3)%d(maxmasses%d(i3))
          maxmasses%d(i3) = maxmasses%d(i3)-1
          partof(j3)%d(j2) = partof(j3)%d(maxmasses%d(j3))
          maxmasses%d(j3) = maxmasses%d(j3)-1
        end if
        if (k1 > coef**k3 - 1) then
          k3 = k3 + 1
        end if
        maxmasses%d(k3) = maxmasses%d(k3)+1
        partof(k3)%d(maxmasses%d(k3)) = k1
        
        u%d(k3) = max(u%d(k3),u_vec(k1))
        v%d(k3) = max(v%d(k3),v_vec(k1))
        
        totali = totali - ratesi%d(i3)
        ratesi%d(i3) = maxmasses%d(i3)*u%d(i3)
        totali = totali + ratesi%d(i3)
        totali = totali - ratesi%d(j3)
        ratesi%d(j3) = maxmasses%d(j3)*u%d(j3)
        totali = totali + ratesi%d(j3)
        totali = totali - ratesi%d(k3)
        ratesi%d(k3) = maxmasses%d(k3)*u%d(k3)
        totali = totali + ratesi%d(k3)
  
        totalj = totalj - ratesj%d(i3)
        ratesj%d(i3) = maxmasses%d(i3)*v%d(i3)
        totalj = totalj + ratesj%d(i3)
        totalj = totalj - ratesj%d(j3)
        ratesj%d(j3) = maxmasses%d(j3)*v%d(j3)
        totalj = totalj + ratesj%d(j3)
        totalj = totalj - ratesj%d(k3)
        ratesj%d(k3) = maxmasses%d(k3)*v%d(k3)
        totalj = totalj + ratesj%d(k3)
      
      totalrate = totali*totalj
      
    end if

    if (np <= maxnp/2) then
      np = np*2
      nvi%d(:) = nvi%d(:)*2
      do j = 1, maxpartof
        do j1 = 1, maxmasses%d(j)
          partof(j)%d(j1+maxmasses%d(j)) = partof(j)%d(j1)
        end do
        maxmasses%d(j) = 2*maxmasses%d(j)
      end do
      ratesi%d(i3min:maxpartof) = 2*ratesi%d(i3min:maxpartof)
      ratesj%d(i3min:maxpartof) = 2*ratesj%d(i3min:maxpartof)
      totali = sum(ratesi%d(i3min:maxpartof))
      totalj = sum(ratesj%d(i3min:maxpartof))
      totalrate = totali*totalj
      nv1 = nv1/2
    end if
    
    if ((verbose > 1) .and. ((mod(step,everys) == 0) .or. (curtime >= maxtime))) then
      print *, step, curtime, np, maxmass
      print *, 'Time passed', dsecnd() - dstart
    end if
  end do
!   close(2)
  nv0 = nv1*np
 
  if (verbose > 0) then
    if (verbose == 1) then
      print *, 'Time passed', dsecnd() - dstart
    end if
    print *, 'Collisions', step
    print *, 'Number density', nv0
    print *, 'Maxmass', maxmass
    print *, 'DONE'
  end if
end

!sigma^2 part of ballistic kernel
subroutine TempApproximate(u, v, ui, vi, nvi, temps, posit)
  Type(Mtrx), intent(out) :: u, v, ui, vi
  Type(IntVec), intent(in) :: nvi
  Type(Vector), intent(in) :: temps
  Integer(4), intent(out) :: posit

  Integer(4) n, i, j, r

  n = nvi%n
  call TempKernel(n, u, v)
  
  r = u%m
  call ui%init(n, r)
  call vi%init(n, r)
  !1/2*C_{ij} = 1/2*sigma_{ij}*(sqrt(T_i/j) + sqrt(T_j/j)) ~ sigma_{ij}/sqrt(i)*sqrt(T_i)
  do i = 1, r
    do j = 1, n
      ui%d(j,i) = nvi%d(j)*u%d(j,i)*Utemp(temps%d(j), j, i)
      vi%d(j,i) = nvi%d(j)*v%d(j,i)*Vtemp(temps%d(j), j, i)
    end do
  end do
  
  !All approximation factors are positive
  posit = 1
end

subroutine TempUpdate(i, u, v, ui, vi, nvec, tvec)
  Integer(4), intent(in) :: i
  Type(Mtrx), intent(in) :: u, v
  Type(Mtrx) :: ui, vi
  Type(IntVec), intent(in) :: nvec
  Type(Vector), intent(in) :: tvec
  
  Integer(4) j

  do j = 1, u%m
    ui%d(i,j) = Utemp(tvec%d(i), i, j)*u%d(i,j)*nvec%d(i)
    vi%d(i,j) = Vtemp(tvec%d(i), i, j)*v%d(i,j)*nvec%d(i)
  end do
end

!Random direction (unit vector)
subroutine randir(ex, ey, ez)
  Double precision, intent(out) :: ex, ey, ez
  Double precision costh, sinth, phi
  phi = grnd()*2*pi
  costh = grnd()*2-1
  sinth = sqrt(1.0d0 - costh**2)
  ex = cos(phi) * sinth
  ey = sin(phi) * sinth
  ez = costh
end

!Random Gaussian vector, mean power 1
subroutine randiradv(ex, ey, ez)
  Double precision, intent(out) :: ex, ey, ez
  call ggrnd2(ex, ey)
  ex = ex/sqrt(3.0d0)
  ey = ey/sqrt(3.0d0)
  ez = ggrnd()/sqrt(3.0d0)
end

function treebuild(mas, treestart) Result(res)
  Type(Mtrx), intent(in) :: mas
  Integer(4), intent(in) :: treestart
  Type(Mtrx) :: res
  Integer(4) n, r, m, i, j
  n = mas%n
  r = mas%m
  m = treestart+ISHFT(n,-2)-1
  call res%init(m,r)
  do j = 1, r
    do i = treestart, m
      res%d(i,j) = mas%d(ISHFT(i-treestart,2)+1,j) + mas%d(ISHFT(i-treestart,2)+2,j) + &
      mas%d(ISHFT(i-treestart,2)+3,j) + mas%d(ISHFT(i-treestart,2)+4,j)
    end do
    do i = treestart-1, 1, -1
      res%d(i,j) = res%d(ISHFT(i,2)+1,j) + res%d(ISHFT(i,2)+2,j) + &
      res%d(ISHFT(i,2)+3,j) + res%d(ISHFT(i,2)+4,j)
    end do
  end do
end

function treebuildv(mas, treestart) Result(res)
  Type(Vector), intent(in) :: mas
  Integer(4), intent(in) :: treestart
  Type(Vector) :: res
  Integer(4) n, m, i
  n = mas%n
  m = treestart+ISHFT(n,-2)-1
  call res%init(m)
  do i = treestart, m
    res%d(i) = mas%d(ISHFT(i-treestart,2)+1) + mas%d(ISHFT(i-treestart,2)+2) + &
    mas%d(ISHFT(i-treestart,2)+3) + mas%d(ISHFT(i-treestart,2)+4)
  end do
  do i = treestart-1, 1, -1
    res%d(i) = res%d(ISHFT(i,2)+1) + res%d(ISHFT(i,2)+2) + &
    res%d(ISHFT(i,2)+3) + res%d(ISHFT(i,2)+4)
  end do
end

subroutine treeupdate2(tree, mas, jin, treestart, upd)
  Type(Mtrx) :: tree
  Type(Mtrx), intent(in) :: mas
  Integer(4), intent(in) :: jin
  Integer(4), intent(in) :: treestart
  Type(Vector), intent(in) :: upd
  Integer(4) n, j
  n = mas%n
  j = treestart + ISHFT(jin-1,-2)
  tree%d(j,:) = tree%d(j,:) + upd%d(:)
  do while (j > 4)
    j = ISHFT(j-1,-2)
    tree%d(j,:) = tree%d(j,:) + upd%d(:)
  end do
end

subroutine treeupdatev(tree, mas, jin, treestart, upd)
  Type(Vector) :: tree
  Type(Vector), intent(in) :: mas
  Integer(4), intent(in) :: jin
  Integer(4), intent(in) :: treestart
  Double precision, intent(in) :: upd
  Integer(4) n, j
  n = mas%n
  j = treestart + ISHFT(jin-1,-2)
  tree%d(j) = tree%d(j) + upd
  do while (j > 4)
    j = ISHFT(j-1,-2)
    tree%d(j) = tree%d(j)+upd
  end do
end

!Find value in tree № curr
function treefindr(tree, mas, randi, curr, treestart) Result(res)
  Type(Mtrx), intent(in) :: tree, mas
  Double precision :: randi
  Integer(4), intent(in) :: curr
  Integer(4), intent(in) :: treestart
  Integer(4) :: res
  Integer(4) n, j, m
  n = mas%n
  m = tree%n
  j = 1
  do while (j < m)
    if (randi < tree%d(j, curr)) then
      j = ISHFT(j,2)+1
      cycle
    end if
    randi = randi - tree%d(j, curr)
    if (randi < tree%d(j+1, curr)) then
      j = ISHFT(j+1,2)+1
      cycle
    end if
    randi = randi - tree%d(j+1, curr)
    if (randi < tree%d(j+2, curr)) then
      j = ISHFT(j+2,2)+1
      cycle
    end if
    randi = randi - tree%d(j+2, curr)
    j = ISHFT(j+3,2)+1
  end do
  j = j-m
  
  if (randi < mas%d(j, curr)) then
    res = j
    return
  end if
  randi = randi - mas%d(j, curr)
  if (randi < mas%d(j+1, curr)) then
    res = j+1
    return
  end if
  randi = randi - mas%d(j+1, curr)
  if (randi < mas%d(j+2, curr)) then
    res = j+2
    return
  end if
  res = j+3
end

function treefindv(tree, mas, randi, treestart) Result(res)
  Type(Vector), intent(in) :: tree, mas
  Double precision :: randi
  Integer(4), intent(in) :: treestart
  Integer(4) :: res
  Integer(4) n, j, m
  n = mas%n
  m = tree%n
  j = 1
  do while (j < m)
    if (randi < tree%d(j)) then
      j = ISHFT(j,2)+1
      cycle
    end if
    randi = randi - tree%d(j)
    if (randi < tree%d(j+1)) then
      j = ISHFT(j+1,2)+1
      cycle
    end if
    randi = randi - tree%d(j+1)
    if (randi < tree%d(j+2)) then
      j = ISHFT(j+2,2)+1
      cycle
    end if
    randi = randi - tree%d(j+2)
    j = ISHFT(j+3,2)+1
  end do
  j = j-m
  
  if (randi < mas%d(j)) then
    res = j
    return
  end if
  randi = randi - mas%d(j)
  if (randi < mas%d(j+1)) then
    res = j+1
    return
  end if
  randi = randi - mas%d(j+1)
  if (randi < mas%d(j+2)) then
    res = j+2
    return
  end if
  res = j+3
end

!Find value in the sum of trees
function treefind(tree, mas, vec, randi, treestart) Result(res)
  Type(Mtrx), intent(in) :: tree, mas
  Type(Vector), intent(in) :: vec
  Double precision :: randi
  Integer(4), intent(in) :: treestart
  Integer(4) :: res
  Double precision val
  Integer(4) n, r, j, m
  n = mas%n
  r = mas%m
  m = tree%n
  j = 1
  do while (j < m)
    val = tovec(tree%subarray(j, r, j, 1))*vec
    if (randi < val) then
      j = ISHFT(j,2)+1
      cycle
    end if
    randi = randi - val
    val = tovec(tree%subarray(j+1, r, j+1, 1))*vec
    if (randi < val) then
      j = ISHFT(j+1,2)+1
      cycle
    end if
    randi = randi - val
    val = tovec(tree%subarray(j+2, r, j+2, 1))*vec
    if (randi < val) then
      j = ISHFT(j+2,2)+1
      cycle
    end if
    randi = randi - val
    j = ISHFT(j+3,2)+1
  end do
  j = j-m
  
  val = tovec(mas%subarray(j, r, j, 1))*vec
  if (randi < val) then
    res = j
    return
  end if
  randi = randi - val
  val = tovec(mas%subarray(j+1, r, j+1, 1))*vec
  if (randi < val) then
    res = j+1
    return
  end if
  randi = randi - val
  if (randi < tovec(mas%subarray(j+1, r, j+1, 1))*vec) then
    res = j+2
    return
  end if
  res = j+3
end

!Integer to string converter
recursive function itoa(i,j) result(res)
  character(:),allocatable :: res
  integer,intent(in) :: i
  integer,intent(in),optional :: j
  character(range(i)+2) :: tmp
  if (present(j)) then
    write(tmp,'(i0.'//itoa(j)//')') i
  else
    write(tmp,'(i0)') i
  end if
  res = trim(tmp)
end function

end module
