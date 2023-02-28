
program bubr
  !Integer(4) mt
  !Integer(4) i
  call init_random_seed()
  call ExampleA()
  call ExampleB()
  call ExampleP()
  call ExampleR()
  !do mt = 1, 4
  !  do i = 0, 20
  !    call reconstruct_test(i*1.0d0, mt)
  !  end do
  !end do
  
  contains

include "incfiles/ExampleA.f90"
include "incfiles/ExampleB.f90"
include "incfiles/ExampleP.f90"
include "incfiles/ExampleR.f90"

!Random seed generator
subroutine init_random_seed()
  integer :: i, n, clock
  integer, dimension(:), allocatable :: seed
  n = 1
  clock = 1

  call random_seed(size = n)
  allocate(seed(n))

  call system_clock(count=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  call random_seed(put = seed)

  deallocate(seed)
end subroutine init_random_seed

!Integer to string converter
function itoa(i) result(res)
  character(:),allocatable :: res
  integer,intent(in) :: i
  character(range(i)+2) :: tmp
  write(tmp,'(i0)') i
  res = trim(tmp)
end function

end program
