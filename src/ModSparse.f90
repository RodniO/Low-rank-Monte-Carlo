Module ModSparse
  Use ModMtrx

  !Sparse matrix type in CSR format
  Type SparseRow
    Integer(4) nz !Number of nonzero elements
    Integer(4) n !Number of rows
    Integer(4), Allocatable :: j(:) !Column numbers for nonzero elements (size nz)
    Integer(4), Allocatable :: i(:) !Total number of nonzero elements till each row (size n+1)
    Real(8), Allocatable :: d(:) !Values of nonzero matrix elements (size nz)
    Contains
      Procedure :: init => Sparse_constructor !Initialize with zeros
      Procedure :: transform => Sparse_transform !Full to sparse transform
      Procedure :: detransform => Sparse_detransform !Sparse to full transform
      Procedure :: deinit => Sparse_destructor !Deallocate all arrays
      Procedure :: mask => Sparse_mask !Returns sparse mask of the sparse matrix
      Procedure :: get => Sparse_get !Get alement from position
      Procedure :: fnorm => Sparse_fnorm
  End type
  
  !Sparse matrix type in CSC format
  Type SparseCol
    Integer(4) nz !Number of nonzero elements
    Integer(4) n !Number of columns
    Integer(4), Allocatable :: i(:) !Row numbers for nonzero elements (size nz)
    Integer(4), Allocatable :: j(:) !Total number of nonzero elements till each column (size n+1)
    Real(8), Allocatable :: d(:) !Values of nonzero matrix elements (size nz)
    Contains
      Procedure :: init => SparseC_constructor !Initialize with zeros
      Procedure :: transform => SparseC_transform !Full to sparse transform
      Procedure :: detransform => SparseC_detransform !Sparse to full transform
      Procedure :: deinit => SparseC_destructor !Deallocate all arrays
  End type
  
  !Create sparse matrix from full matrix
  interface assignment(=)
    module procedure Sparse_transform, SparseC_sparse
  end interface
  
  !Multiply by vector or full matrix
  interface operator(*)
    module procedure Sparse_mvec, Sparse_mmat
  end interface

  Contains
  
    subroutine Sparse_constructor(this, nz, n)
      Class(SparseRow) :: this
      Integer(4) nz
      Integer(4) n, i
      this%nz = nz
      this%n = n
      Allocate(this%d(nz))
      Allocate(this%j(nz))
      Allocate(this%i(n+1))
      do i = 1, nz
        this%d(i) = 0
        this%j(i) = 0
      end do
      do i = 1, n+1
        this%i(i) = 0
      end do
    end
    
    subroutine SparseC_constructor(this, nz, n)
      Class(SparseCol) :: this
      Integer(4) nz
      Integer(4) n, i
      this%nz = nz
      this%n = n
      Allocate(this%d(nz))
      Allocate(this%i(nz))
      Allocate(this%j(n+1))
      do i = 1, nz
        this%d(i) = 0
        this%i(i) = 0
      end do
      do i = 1, n+1
        this%j(i) = 0
      end do
    end
    
    subroutine Sparse_transform(this, mat)
      Type(Mtrx), intent(in) :: mat
      Class(SparseRow), intent(out) :: this
      Integer(4) nz, n, m
      Integer(4) i, j
      n = mat%n
      m = mat%m
      this%n = n
      Allocate(this%i(n+1))
      nz = 0
      do i = 1,n
        this%i(i) = nz+1
        do j = 1,m
          if (mat%d(i, j) /= 0.0d0) then
            nz = nz + 1
          endif
        end do
      end do
      this%i(n+1) = nz+1
      this%nz = nz
      Allocate(this%d(nz))
      Allocate(this%j(nz))
      nz = 0
      do i = 1,n
        do j = 1,m
          if (mat%d(i, j) /= 0.0d0) then
            nz = nz + 1
            this%d(nz) = mat%d(i, j)
            this%j(nz) = j
          endif
        end do
      end do
    end
    
    subroutine SparseC_transform(this, mat)
      Type(Mtrx), intent(in) :: mat
      Class(SparseCol), intent(out) :: this
      Integer(4) nz, n, m
      Integer(4) i, j
      n = mat%n
      m = mat%m
      this%n = m
      Allocate(this%j(m+1))
      nz = 0
      do j = 1,m
        this%j(j) = nz+1
        do i = 1,n
          if (mat%d(i, j) /= 0.0d0) then
            nz = nz + 1
          endif
        end do
      end do
      this%j(n+1) = nz+1
      this%nz = nz
      Allocate(this%d(nz))
      Allocate(this%i(nz))
      nz = 0
      do j = 1,m
        do i = 1,n
          if (mat%d(i, j) /= 0.0d0) then
            nz = nz + 1
            this%d(nz) = mat%d(i, j)
            this%i(nz) = i
          endif
        end do
      end do
    end
    
    subroutine SparseC_sparse(this, mat)
      Type(SparseRow), intent(in) :: mat
      Class(SparseCol), intent(out) :: this
      Integer(4) nz, n, m
      Integer(4) i, j
      Integer(4), allocatable :: counters(:)
      n = mat%n
      m = maxval(mat%j)
      this%n = m
      Allocate(this%j(m+1))
      nz = mat%nz
      this%nz = nz
      Allocate(this%d(nz))
      Allocate(this%i(nz))
      Allocate(counters(m))
      
      do i = 1, m+1
        this%j(i) = 0
      end do
      do i = 1, mat%nz
        this%j(mat%j(i)+1) = this%j(mat%j(i)+1) + 1
      end do
      do i = 1, m
        this%j(i+1) = this%j(i) + this%j(i+1)
        this%j(i) = this%j(i)+1
        counters(i) = this%j(i)
      end do
      this%j(m+1) = nz+1

      j = 0
      do nz = 1, mat%nz
        do while (mat%i(j+1) <= nz)
          j = j + 1
        end do
        this%d(counters(mat%j(nz))) = mat%d(nz)
        this%i(counters(mat%j(nz))) = j
        counters(mat%j(nz)) = counters(mat%j(nz)) + 1
      end do
    end
    
    function Sparse_detransform(this, n, perm) Result(res)
      Class(SparseRow) :: this
      Integer(4), intent(in), optional :: n !rows
      Type(IntVec), optional :: perm
      Type(Mtrx) :: res
      Integer(4) m, n_, i, j, nz
      m = 0
      j = 0
      do i = 1,this%nz
        m = max(m, this%j(i))
      end do
      if (present(n)) then
        n_ = n
      else
        n_ = this%n
      end if
      call res%init(n_, m)
      if (present(perm)) then
        do i = 1, n_
          j = perm%d(i)
          do nz = this%i(j), this%i(j+1)-1
            res%d(i, this%j(nz)) = this%d(nz)
          end do
        end do
      else
        do nz = 1, this%i(n_+1)-1
          do while (this%i(j+1) <= nz)
            j = j + 1
          end do
          res%d(j, this%j(nz)) = this%d(nz)
        end do
      end if
    end
    
    function SparseC_detransform(this, m, perm) Result(res)
      Class(SparseCol) :: this
      Integer(4), intent(in), optional :: m !columns
      Type(IntVec), optional :: perm
      Type(Mtrx) :: res
      Integer(4) m_, n, i, j, nz
      n = 0
      j = 0
      do i = 1,this%nz
        n = max(n, this%i(i))
      end do
      if (present(m)) then
        m_ = m
      else
        m_ = this%n
      end if
      call res%init(n, m_)
      if (present(perm)) then
        do i = 1, m_
          j = perm%d(i)
          do nz = this%j(j), this%j(j+1)-1
            res%d(this%i(nz), i) = this%d(nz)
          end do
        end do
      else
        do nz = 1, this%j(m_+1)-1
          do while (this%j(j+1) <= nz)
            j = j + 1
          end do
          res%d(this%i(nz), j) = this%d(nz)
        end do
      end if
    end
    
    function Sparse_mvec(this, v) Result(res)
      Type(SparseRow), intent(in) :: this
      Type(Vector), intent(in) :: v
      Type(vector) :: res
      Integer(4) j, nz
      call res%init(this%n)
      j = 0
      do nz = 1, this%nz
        do while (this%i(j+1) <= nz)
          j = j + 1
        end do
        res%d(j) = res%d(j) + this%d(nz) * v%d(this%j(nz))
      end do
    end
    
    function Sparse_mmat(this, mat) Result(res)
      Type(SparseRow), intent(in) :: this
      Type(Mtrx), intent(in) :: mat
      Type(Mtrx) :: res
      Integer(4) i, j, nz
      call res%init(this%n, mat%m)
      do i = 1, mat%m
        j = 0
        do nz = 1, this%nz
          do while (this%i(j+1) <= nz)
            j = j + 1
          end do
          res%d(j, i) = res%d(j, i) + this%d(nz) * mat%d(this%j(nz), i)
        end do
      end do
    end
    
    function Sparse_mask(this) Result(res)
      Class(SparseRow), intent(in) :: this
      Type(SparseRow) :: res
      
      call res%init(this%nz, this%n)
      res%d(:) = 1
      res%i(:) = this%i(:)
      res%j(:) = this%j(:)
    end
    
    function Sparse_get(this, i, j) Result(res)
      Class(SparseRow), intent(in) :: this
      Integer(4), intent(in) :: i, j
      Double precision :: res
      Integer(4) k
      
      res = 0
      do k = this%i(i), this%i(i+1)-1
        if (this%j(k) == j) then
          res = this%d(k)
        end if
      end do
    end
    
    function Sparse_fnorm(this) Result(res)
      Class(SparseRow), intent(in) :: this
      Double precision :: res
      Double precision dnrm2
      res = dnrm2(this%nz, this%d, 1)
    end
    
    subroutine Sparse_destructor(this)
      Class(SparseRow) :: this
      Deallocate(this%d)
      Deallocate(this%i)
      Deallocate(this%j)
    end
    
    subroutine SparseC_destructor(this)
      Class(SparseCol) :: this
      Deallocate(this%d)
      Deallocate(this%i)
      Deallocate(this%j)
    end
end
