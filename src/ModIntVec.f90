Module ModIntVec

  !Module for vector of integers
  
  !Useful constants
  DOUBLE PRECISION, parameter :: eps = 1.0d-15
  DOUBLE PRECISION, parameter :: pi = acos(-1.0d0)
  
  !Vector type
  Type IntVec
    Integer(4) n !Dimension
    Integer(4),Allocatable :: d(:) !Array of vector elements of size n
    Contains
      Procedure :: init => intvec_constructor !Allocates and initializes with zeros
      Procedure :: deinit => intvec_destructor !Deallocates
      Procedure :: set => intvec_set !Sets d to desired array
      Procedure :: fromsubset => intvec_fromsubset !Creates permutation from its subarray
      Procedure :: extend => intvec_extend !Extends subarray to permutation
      Procedure :: permapp => intvec_permapp !Apply permutation
      Procedure :: perm => intvec_perm !Initialize identity permutation
      Procedure :: subarray => intvec_subarray !Returns subvector
      Procedure :: bst => intvec_bst !Binary search tree order
      Procedure :: cheb => intvec_cheb !Chebyshev polynomial zeros order
      Generic :: permrand => intvec_permrand, intvec_permrandk !Random permutation
      Procedure, private :: intvec_permrand, intvec_permrandk
      Procedure :: swap => intvec_swap !Swap two vector elements
      Procedure :: reverse => intvec_reverse !Reverse order of elements
      Procedure :: copy => intvec_copy !Copy vector
      Procedure :: perminv => intvec_perminv !Inverse permutation
      Procedure :: sort => ivec_sort !Vector sort
  End type
  
  !Turns integer array to IntVec type
  interface assignment(=)
    module procedure IntVecArray_transform
  end interface
  
  Contains
  !See brief descriptions above

    subroutine ivec_sort(this, per)
      Class(IntVec) :: this
      Type(IntVec) :: per
      call iintvec_sort(this%n, per%d, this%d)
    end
    
    function intvec_subarray(this, n2, n1) Result(res)
      Class(IntVec), intent(in) :: this
      Integer(4), intent(in) :: n2 !End index
      Integer(4), intent(in), optional :: n1 !Start index
      Integer(4) n1_
      Type(IntVec) res
      if (present(n1)) then
        n1_ = n1
      else
        n1_ = 1
      end if
      res%n = n2 - n1_ + 1
      Allocate(res%d(res%n))
      res%d(:) = this%d(n1_:n2)
    end

    subroutine intvec_permapp(this, per, c)
      Class(IntVec) :: this
      Type(IntVec), intent(in) :: per
      Integer(4), intent(in) :: c
      Type(IntVec) :: res
      Integer(4) i
      call res%copy(this)
      if (c == 1) then
        do i = 1, per%n
          this%d(i) = res%d(per%d(i))
        end do
      else
        do i = 1, per%n
          this%d(per%d(i)) = res%d(i)
        end do
      end if
    end
    
    subroutine intvec_swap(this, a, b)
      Class(IntVec) :: this
      Integer(4), intent(in) :: a, b
      Integer(4) tmp
      !if ((a > this%n) .or. (b > this%n)) then
      !  print *, "error in swap_vec", a, b, this%n
        !call backtrace()
      !  return
      !end if
      tmp = this%d(a)
      this%d(a) = this%d(b)
      this%d(b) = tmp
    end
    
    subroutine intvec_permrand(this, n)
      Class(IntVec) :: this
      Integer(4), intent(in) :: n
      
      Double precision :: vec(n)
      Integer(4) i, j
      
      call this%perm(n)
      call random_number(vec)
      do i = 1, n-1
        j = i + floor((n-i+1)*vec(i))
        call this%swap(i,j)
      end do
    end
    
    subroutine intvec_permrandk(this, n, k)
      Class(IntVec) :: this
      Integer(4), intent(in) :: n
      Integer(4), intent(in) :: k
      
      Double precision :: vec(k)
      Integer(4) i, j
      
      call this%perm(n)
      call random_number(vec)
      do i = 1, k
        j = i + floor((n-i+1)*vec(i))
        call this%swap(i,j)
      end do
    end
    
    function intvec_reverse(this) Result(res)
      Class(IntVec) :: this
      Type(IntVec) :: res
      Integer(4) i
      res%n = this%n
      Allocate(res%d(res%n))
      do i = 1, res%n
        res%d(i) = this%d(this%n-i+1)
      end do
    end
    
    subroutine IntVecArray_transform(this, array)
      Integer(4), dimension(:), intent(in) :: array
      Class(IntVec), intent(out) :: this
      this%n = size(array)
      Allocate(this%d(this%n))
      this%d(:) = array(:)
    end
    
    subroutine intvec_copy(this, v)
      Class(IntVec) :: this
      Type(IntVec), intent(in) :: v
      if (.not. allocated(this%d)) then
        Allocate(this%d(v%n))
      else if (this%n < v%n) then
        Deallocate(this%d)
        Allocate(this%d(v%n))
      end if
      this%n = v%n
      this%d(1:this%n) = v%d(1:this%n)
    end
    
    subroutine intvec_perm(this, n)
      Class(IntVec) :: this
      Integer(4), intent(in) :: n
      Integer(4) i
      this%n = n
      if (allocated(this%d)) then
        Deallocate(this%d)
      end if
      Allocate(this%d(n))
      do i = 1, n
        this%d(i) = i
      end do
    end
    
    subroutine intvec_bst(this, n)
      Class(IntVec) :: this
      Integer(4), intent(in) :: n
      Integer(4), allocatable :: lengths(:)
      Integer(4) i, ind
      this%n = n
      if (allocated(this%d)) then
        Deallocate(this%d)
      end if
      Allocate(this%d(n))
      Allocate(lengths(n))
      ind = 1
      i = 1
      this%d(1) = (n+1)/2
      lengths(1) = n
      do while (ind < n)
        if (lengths(i) > 2) then
          ind = ind + 1
          this%d(ind) = this%d(i) - (lengths(i)+3)/4
          lengths(ind) = (lengths(i)-1)/2
        end if
        if (lengths(i) > 1) then
          ind = ind + 1
          this%d(ind) = this%d(i) + (lengths(i)+2)/4
          lengths(ind) = lengths(i)/2
        end if
        i = i + 1
      end do
    end
    
    subroutine intvec_cheb(this, n)
      Class(IntVec) :: this
      Integer(4), intent(in) :: n
      Type(IntVec) bst
      Double precision, allocatable :: cheb(:)
      logical, allocatable :: used(:)
      Double precision phi
      Integer(4) i, j, ind, nzeros
      
      this%n = n
      if (allocated(this%d)) then
        Deallocate(this%d)
      end if
      Allocate(this%d(n))
      Allocate(used(n))
      used = .false.
      nzeros = floor((pi/4)*n)
      Allocate(cheb(nzeros))
      phi = pi/(nzeros+1)
      do i = 1, nzeros
        cheb(i) = (1-cos(i*phi))/2
      end do
      call bst%bst(nzeros)
      ind = 1
      do i = 1, nzeros
        j = floor(n*cheb(bst%d(i)))+1
        if (.not.used(j)) then
          this%d(ind) = j
          used(j) = .true.
          ind = ind + 1
        end if
      end do
      do i = 1, n
        if (.not.used(i)) then
          this%d(ind) = i
          ind = ind + 1
        end if
      end do
    end
  
    subroutine intvec_constructor(this, n)
      Class(IntVec) :: this
      Integer(4) n
      this%n = n
      if (allocated(this%d)) then
        Deallocate(this%d)
      end if
      Allocate(this%d(n))
      this%d(:) = 0
    end
    
    subroutine intvec_set(this, d)
      Class(IntVec) :: this
      Integer(4) :: d(:)
      this%d = d
      this%n = size(d)
    end
    
    subroutine intvec_destructor(this)
      Class(IntVec) :: this
      this%n = 0
      Deallocate(this%d)
    end
    
    function intvec_perminv(this) Result(res)
      Class(IntVec) :: this
      Type(IntVec) :: res
      Integer(4) i
      res%n = this%n
      Allocate(res%d(this%n))
      do i = 1, this%n
        res%d(this%d(i)) = i
      end do
    end
    
    subroutine iintvec_sort(array_size,index,value)
      integer, intent(in) :: array_size
      integer, intent(inout) :: index(array_size)
      Integer(4), intent(in) :: value(array_size)
      Integer(4) :: QSORT_THRESHOLD = 8
      include "qsort_inline.inc"
    contains
      include "qsort_inline_index.inc"
      logical &
      function less_than(a,b)
        integer, intent(in) :: a,b
        if ( value(index(a)) == value(index(b)) ) then
          less_than = index(a) < index(b)
        else
          less_than = value(index(a)) < value(index(b))
        end if
      end function
    end subroutine
    
    subroutine intvec_fromsubset(this, n, k, d)
      Class(IntVec) :: this
      Integer(4), intent(in) :: n, k
      Integer(4), intent(in) :: d(k)

      Integer(4) i
      Integer(4) swappedto(k)
      
      call this%perm(n)
      do i = 1, k
        swappedto(i) = i
      end do
      do i = 1, k
        if (d(i) > k) then
          swappedto(this%d(i)) = d(i)
          call this%swap(i,d(i))
        else
          swappedto(this%d(i)) = swappedto(d(i))
          call this%swap(i,swappedto(d(i)))
          swappedto(d(i)) = i
        end if
      end do
    end
    
    subroutine intvec_extend(this, n)
      Class(IntVec) :: this
      Integer(4), intent(in) :: n

      Integer(4), allocatable :: newthis(:)
      Integer(4) i, k, tmp
      Integer(4), allocatable :: swappedto(:)
      k = this%n
      this%n = n
      Allocate(newthis(n))
      do i = 1, n
        newthis(i) = i
      end do
      Allocate(swappedto(k))
      do i = 1, k
        swappedto(i) = i
      end do
      do i = 1, k
        tmp = newthis(i)
        if (this%d(i) > k) then
          newthis(i) = newthis(this%d(i))
          newthis(this%d(i)) = tmp
          swappedto(tmp) = this%d(i)
        else
          newthis(i) = newthis(swappedto(this%d(i)))
          newthis(swappedto(this%d(i))) = tmp
          swappedto(tmp) = swappedto(this%d(i))
          swappedto(this%d(i)) = i
        end if
      end do
      this%d = newthis
    end
    
    !Create vector of ones or k-th vector of standard basis
    elemental function eveci(n, k) Result(res)
      Integer(4), intent(in) :: n
      Integer(4), intent(in), optional :: k
      Type(IntVec) :: res
      res%n = n
      Allocate(res%d(n))
      if (present(k)) then
        res%d(:) = 0
        res%d(k) = 1
      else
        res%d(:) = 1
      end if
    end

end
