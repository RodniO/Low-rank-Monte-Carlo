Module ModIntVec

  !Module for vector of integers
  
  !Vector type
  Type IntVec
    Integer(4) n !Dimension
    Integer(4),Allocatable :: d(:) !Array of vector elements of size n
    Contains
      Procedure :: init => intvec_constructor !Allocates and initializes with zeros
      Procedure :: deinit => intvec_destructor !Deallocates
      Procedure :: set => intvec_set !Sets d to desired array
      Procedure :: permapp => intvec_permapp !Apply permutation
      Procedure :: perm => intvec_perm !Initialize identity permutation
      Generic :: permrand => intvec_permrand, intvec_permrandk !Random permutation
      Procedure, private :: intvec_permrand, intvec_permrandk
      Procedure :: swap => intvec_swap !Swap two vector elements
      Procedure :: reverse => intvec_reverse !Reverse order of elements
      Procedure :: copy => intvec_copy !Copy vector
      Procedure :: perminv => intvec_perminv !Inverse permutation
  End type
  
  !Turns integer array to IntVec type
  interface assignment(=)
    module procedure IntVecArray_transform
  end interface
  
  Contains
  !See brief descriptions above

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
      if ((a > this%n) .or. (b > this%n)) then
        print *, "error in swap_vec", a, b, this%n
        !call backtrace()
        return
      end if
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
  
    subroutine intvec_constructor(this, n)
      Class(IntVec) :: this
      Integer(4) n
      this%n = n
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
