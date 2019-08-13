
program bubr
  call init_random_seed()
  call ExampleA()
  call ExampleB()
  
  contains

include "incfiles/ExampleA.f90"
include "incfiles/ExampleB.f90"

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
