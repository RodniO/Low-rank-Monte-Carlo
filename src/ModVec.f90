Module ModVec
  USE ModIntVec

  !Wrapper for BLAS vector subroutines
  
  !Vector type
  Type Vector
    Integer(4) n !Dimension
    DOUBLE PRECISION,Allocatable :: d(:) !Array of vector elements of size n
    Contains
      Procedure :: init => Vector_constructor !Allocates and initializes with zeros
      Procedure :: deinit => Vector_destructor !Deallocates
      Procedure :: norm => vec_norm !Euclidean vector norm
      Procedure :: cnorm => vec_cnorm !Infinity (max) norm
      Procedure :: cnormloc => vec_cnormloc !Infinity (max) norm location
      Procedure :: random => Vec_random !Random unit vector
      Procedure :: badrandom => Vec_badrandom !Vector with uniform random elements
      Procedure :: subarray => Vec_subarray !Returns subvector
      Procedure :: permapp => Vec_permapp !Apply permutation
      Procedure :: swap => Vec_swap !Swap two vector elements
      Procedure :: normalize => Vec_normalize !Normalize vector to unity
      Procedure :: maxelement => Vec_maxelement !Find position of maximum element
      Procedure :: house => Vec_house !Create Householder reflection
      Procedure :: sort => Vec_sort !Vector sort
      Procedure :: intarray => Vec_toint !Round values to integers
      Procedure :: reverse => vec_reverse !Reverse order of elements
      Procedure :: shift => vec_shift !Shift elements and put nulls
      Procedure :: copy => vec_copy !Copy vector
      Procedure :: update1 => vec_update1 !Update vector x := ax + y
  End type
  
  !Elementwise product
  interface operator (.dot.)
    module procedure vec_dotmul
  end interface
  
  !Elementwise division
  interface operator (.dd.)
    module procedure vec_dotdiv
  end interface
  
  !Turns integer or double precision array to vector type
  interface assignment(=)
    module procedure IntArray_transform, VecArray_transform
  end interface
  
  !Sum of vectors
  interface operator(+)
    module procedure vec_sum
  end interface
  
  !Vector subtraction
  interface operator(-)
    module procedure vec_sub
  end interface
  
  !Scalar product and vector-number products
  interface operator(*)
    module procedure mul_vec, mul_tnum, mul_tnumr, mul_numt, mul_numtr
  end interface
  
  !Vector-number division
  interface operator(/)
    module procedure div_num, div_numr
  end interface
  
  Contains
  !See brief descriptions above
  
    subroutine vec_sort(this, per)
      Class(Vector) :: this
      Type(IntVec) :: per
      call intvec_sort(this%n, per%d, this%d)
    end
    
    subroutine vec_permapp(this, per, c)
      Class(Vector) :: this
      Type(IntVec), intent(in) :: per
      Integer(4), intent(in) :: c
      Type(Vector) :: res
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
    
    subroutine vec_swap(this, a, b)
      Class(Vector) :: this
      Integer(4), intent(in) :: a, b
      DOUBLE PRECISION tmp
      tmp = this%d(a)
      this%d(a) = this%d(b)
      this%d(b) = tmp
    end
    
    function vec_shift(this, sh, start) Result(res)
      Class(Vector) :: this
      Integer(4), intent(in) :: sh !Shift size (can be negative)
      Integer(4), optional :: start !Shift starts at this index [default 1]
      Type(Vector) :: res !Shifted vector, with zeros in the beginning (for sh > 0)
      Integer(4) start_
      Integer(4) i
      start_ = 1
      if (present(start)) then
        start_ = start
      end if
      call res%init(this%n)
      if (sh >= 0) then
        call dcopy(start_-1, this%d, 1, res%d, 1)
        do i = this%n, start_+sh, -1
          res%d(i) = this%d(i-sh)
        end do
      else
        do i = this%n+sh-start_+2, this%n
          res%d(i) = this%d(i)
        end do
        do i = 1, this%n+sh-start_+1
          res%d(i) = this%d(i-sh)
        end do
      end if
    end
    
    function vec_reverse(this) Result(res)
      Class(Vector) :: this
      Type(Vector) :: res
      Integer(4) i
      res%n = this%n
      Allocate(res%d(res%n))
      do i = 1, res%n
        res%d(i) = this%d(this%n-i+1)
      end do
    end
    
    subroutine IntArray_transform(this, array)
      Integer(4), dimension(:), intent(in) :: array
      Class(Vector), intent(out) :: this
      Integer(4) i
      this%n = size(array)
      Allocate(this%d(this%n))
      do i = 1, this%n
        this%d(i) = array(i)
      end do
    end
    
    subroutine VecArray_transform(this, array)
      DOUBLE PRECISION, dimension(:), intent(in) :: array
      Class(Vector), intent(out) :: this
      this%n = size(array)
      this%d = array
    end
    
    subroutine vec_copy(this, v)
      Class(Vector) :: this
      Type(Vector), intent(in) :: v
      if (.not. allocated(this%d)) then
        Allocate(this%d(v%n))
      else if (this%n < v%n) then
        Deallocate(this%d)
        Allocate(this%d(v%n))
      end if
      this%n = v%n
      call dcopy(v%n, v%d, 1, this%d, 1)
    end
    
    subroutine vec_update1(this, alpha, x)
      Class(Vector) :: this
      DOUBLE PRECISION, intent(in) :: alpha
      Type(Vector), intent(in) :: x
      call daxpy(min(this%n,x%n), alpha, x%d, 1, this%d, 1)
    end
    
    function Vec_toint(this) Result(array)
      Class(Vector), intent(in) :: this
      Integer(4) :: array(this%n)
      Integer(4) i
      do i = 1, this%n
        array(i) = floor(this%d(i) + 0.5)
      end do
    end
  
    subroutine Vec_normalize(this)
      Class(Vector) :: this
      call dscal(this%n, 1.0d0/this%norm(), this%d, 1)
    end
    
    function Vec_cnormloc(this) Result(k)
      Class(Vector) :: this
      Integer(4) k
      k = myidamax(this%n,this%d)
    end
    
    function Vec_maxelement(this) Result(k)
      Class(Vector) :: this
      Integer(4) k
      k = myidmax(this%n,this%d)
    end
    
    subroutine Vec_house(this, v, beta)
      Class(Vector) :: this
      Type(Vector), intent(out) :: v
      DOUBLE PRECISION, intent(out) :: beta
      DOUBLE PRECISION sigma
      sigma = this%norm()
      sigma = sigma**2 - this%d(1)**2
      call v%copy(this)
      v%d(1) = 1.0d0
      if ((sigma == 0) .and. (this%d(1) >= 0)) then
        beta = 0.0d0
      else if ((sigma == 0) .and. (this%d(1) < 0)) then
        beta = -2.0d0
      else
        v%d(1) = this%d(1) - this%norm()
        beta = 2.0d0*v%d(1)**2/(sigma + v%d(1)**2)
        v = v/v%d(1)
      end if
    end
    
    subroutine Vec_random(this, n, nz)
      Class(Vector) :: this
      Integer(4) n, i
      Integer(4), Optional, intent(in) :: nz
      DOUBLE PRECISION r1, r2, r
      this%n = n
      Allocate(this%d(n))
      do i = 1, n
        r = 1
        do while (r >= 1)
          call random_number(r1)
          call random_number(r2)
          r1 = r1 * 2 - 1
          r2 = r2 * 2 - 1
          r = r1 * r1 + r2 * r2
        end do
        r2 = r1 * sqrt((-2) * log(r) / r)
        this%d(i) = r2
      end do
      if (.not. (present(nz))) then
        call this%normalize()
      end if
    end
    
    subroutine Vec_badrandom(this, n)
      Class(Vector) :: this
      Integer(4) n, i
      DOUBLE PRECISION r
      this%n = n
      Allocate(this%d(n))
      do i = 1, n
        call random_number(r)
        this%d(i) = 2*r-1
      end do
    end
    
    !Would be good to enable permutations as input
    !Or create another subroutine (same name?) for that
    function Vec_subarray(this, n2, n1) Result(res)
      Class(Vector), intent(in) :: this
      Integer(4), intent(in) :: n2 !End index
      Integer(4), intent(in), optional :: n1 !Start index
      Integer(4) n1_
      Type(Vector) res
      if (present(n1)) then
        n1_ = n1
      else
        n1_ = 1
      end if
      res%n = n2 - n1_ + 1
      Allocate(res%d(res%n))
      res%d(:) = this%d(n1_:n2)
    end
  
    subroutine Vector_constructor(this, n)
      Class(Vector) :: this
      Integer(4) n
      this%n = n
      if (allocated(this%d)) then
        Deallocate(this%d)
      end if
      Allocate(this%d(n))
      this%d(:) = 0
    end
    
    subroutine Vector_destructor(this)
      Class(Vector) :: this
      this%n = 0
      Deallocate(this%d)
    end
    
    function mul_vec(this, v2) Result(res)
      Type(Vector), intent(in) :: this
      DOUBLE PRECISION res, ddot
      Type(Vector), intent(in) :: v2
      !if (this%n == v2%n) then
        res = ddot(this%n, this%d, 1, v2%d, 1)
      !else
      !  print *, "error mult_vec"
      !endif
    end
      
    elemental function mul_tnum(this, num) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      DOUBLE PRECISION, intent(in) :: num
      res%n = this%n
      res%d = num * this%d
    end
    
    elemental function mul_tnumr(this, num) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      Real(4), intent(in) :: num
      res%n = this%n
      res%d = num * this%d
    end
    
    elemental function mul_numt(num, this) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      DOUBLE PRECISION, intent(in) :: num
      res%n = this%n
      res%d = num * this%d
    end
    
    elemental function mul_numtr(num, this) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      Real(4), intent(in) :: num
      res%n = this%n
      res%d = num * this%d
    end
    
    elemental function div_num(this, num) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      DOUBLE PRECISION, intent(in) :: num
      res = this * (1.0d0 / num)
    end
    
    elemental function div_numr(this, num) Result(res)
      Type(Vector), intent(in) :: this
      Type(Vector) :: res
      Real(4), intent(in) :: num
      res = this * (1.0d0 / num)
    end
    
    function vec_norm(this) Result(res)
      Class(Vector) :: this
      DOUBLE PRECISION res, dnrm2
      res = dnrm2(this%n, this%d, 1)
    end
    
    function vec_cnorm(this) Result(res)
      Class(Vector) :: this
      DOUBLE PRECISION res
      res = abs(this%d(myidamax(this%n, this%d)))
    end
    
    elemental function vec_sum(v1, v2) Result(res)
      Type(Vector), intent(in) :: v1, v2
      Type(Vector) :: res
      res%n = v1%n
      !if (v1%n == v2%n) then
        res%d = v1%d + v2%d
      !else
      !  print *, "error sum_vec", v1%n, v2%n
      !endif
    end
    
    elemental function vec_sub(v1, v2) Result(res)
      Type(Vector), intent(in) :: v1, v2
      Type(Vector) :: res
      res%n = v1%n
      !if (v1%n == v2%n) then
        res%d = v1%d - v2%d
      !else
      !  print *, "error sub_vec", v1%n, v2%n
      !endif
    end
    
    elemental function vec_dotmul(v1, v2) Result(res)
      Type(Vector), intent(in) :: v1, v2
      Type(Vector) :: res
      res%n = v1%n
      !if (v1%n == v2%n) then
        res%d = v1%d * v2%d
      !else
      !  print *, "error dotmul_vec", v1%n, v2%n
      !endif
    end
    
    elemental function vec_dotdiv(v1, v2) Result(res)
      Type(Vector), intent(in) :: v1, v2
      Type(Vector) :: res
      res%n = v1%n
      !if (v1%n == v2%n) then
        res%d = v1%d / v2%d
      !else
      !  print *, "error dotdiv_vec", v1%n, v2%n
      !endif
    end
    
    subroutine intvec_sort(array_size,index,value)
      integer, intent(in) :: array_size
      integer, intent(inout) :: index(array_size)
      DOUBLE PRECISION, intent(in) :: value(array_size)
      Integer(4) :: QSORT_THRESHOLD = 8
      include "qsort_inline.inc"
    contains
      include "qsort_inline_index.inc"
      logical &
      function less_than(a,b)
        integer, intent(in) :: a,b
        real(4), parameter :: small=1.0e-6
        if ( abs(value(index(a))-value(index(b))) < small ) then
          less_than = index(a) < index(b)
        else
          less_than = value(index(a)) < value(index(b))
        end if
      end function
    end subroutine
    
    function myidamax(n, A) Result(res)
      Integer(4) n
      Double precision, intent(in) :: A(n)
      Integer(4) res
  
      Double precision M, x, ab(16)
      Integer(4) k, kf, ij(1)
  
      M = 0
      do k = 1, n-16, 16
        ab(:) = abs(A(k:k+15))
        x = maxval(ab)
        if (x > M) then
          ij = maxloc(ab)+k-1
          M = x
        end if
      end do
      kf = n - k + 1
      ab(1:kf) = abs(A(k:n))
      x = maxval(ab(1:kf))
      if (x > M) then
        ij = maxloc(ab(1:kf))+k-1
        M = x
      end if
      res = ij(1)
    end function
    
    function myidmax(n, A) Result(res)
      Integer(4) n
      Double precision, intent(in) :: A(n)
      Integer(4) res
  
      Double precision M, x, ab(16)
      Integer(4) k, kf, ij(1)
  
      ij(1) = 1
      M = A(1)
      do k = 1, n-16, 16
        ab(:) = A(k:k+15)
        x = maxval(ab)
        if (x > M) then
          ij = maxloc(ab)+k-1
          M = x
        end if
      end do
      kf = n - k + 1
      ab(1:kf) = A(k:n)
      x = maxval(ab(1:kf))
      if (x > M) then
        ij = maxloc(ab(1:kf))+k-1
        M = x
      end if
      res = ij(1)
    end function
    
    !Create vector of ones or k-th vector of standard basis
    elemental function evec(n, k) Result(res)
      Integer(4), intent(in) :: n
      Integer(4), intent(in), optional :: k
      Type(Vector) :: res
      res%n = n
      Allocate(res%d(n))
      if (present(k)) then
        res%d(:) = 0
        res%d(k) = 1.0d0
      else
        res%d(:) = 1.0d0
      end if
    end
    
end
