
!Matrix reconstruction example.
subroutine ExampleR()
  USE ModRecon
  Integer(4) n, k !Size and rank
  Double precision q !Sparsity
  Double precision time !Evaluation time
  Type(Mtrx) mat !Matrix to reconstruct
  Type(Mtrx) U, V !SVD factors
  Type(Vector) S !SVD diagonal
  Type(Mtrx) mask !Sparse mask (zeros and ones)
  Type(SparseRow) KnownSparse !Sparse matrix of known elements
  Type(Mtrx) matmask !Full version of KnownSparse
  Type(Mtrx) E !Error
  Integer(4) i !Counter
  
  !Size
  n = 1000
  !Rank
  k = 10
  !Sparsity
  q = 0.2d0
  !Time
  time = 5.0d0
  
  !WELCOMING AND INITIALIZATION
  print *, ''
  print *, 'Finally, let us perform Matrix Completion!'
  print *, ''
  write(*,'(A,I0,A)') ' We are going to try reconstructing random rank ', k, ' matrix'
  print *, 'with the knowledge of only a part of its elements.'
  print *, ''
  
  !Generating random rank k matrix
  call U%random(n,k)
  call V%random(n,k)
  call S%init(k)
  do i = 1, k
    S%d(i) = 2.0d0/2.0d0**i
  end do
  mat = (U .dot. S)*(.T.V)
  call U%deinit()
  call S%deinit()
  call V%deinit()
  
  write(*,'(A,G0.14)') ' Let us say that we know an element with probability ', q
  
  !Leave elements nonzero with probability q
  call mask%mask(n, q)
  matmask = mask .dot. mat
  KnownSparse = matmask
  
  !Initialize low-rank factors
  call U%init(n,1)
  call V%init(1,n)
  
  write(*,'(A,G0.14,A)') ' and give ourselves ', time, ' seconds...'
  print *, ''
  
  call reconstruct(KnownSparse, k, U, V, time, eps, 2)
  
  print *, ''
  print *, 'We arrive at the relative Frobenius norm error equal to'
  E = mat - U*V
  print *, E%fnorm()/mat%fnorm()
  
  print *, ''
  print *, 'All done.'
  print *, 'Goodbye!'
end

!Test subroutine used to get the data from the paper (in preparation) and
!presented on the 5-th International Conference on Matrix Methods in
!Mathematics and applications (MMMA-2019) by Sergey Petrov:
!"Singular value Projection for Matrix Complition with fast approximate projectors"
subroutine reconstruct_test(alpha, mt)
  USE ModRecon
  Double precision, intent(in) :: alpha !Parameter, affecting sparsity (q)
  !alpha = 0..20 in paper
  Integer(4), intent(in) :: mt !Matrix type
  !Rank of all matrices is r (see below)
  !-1 - multiplication of two matrices with uniformly distributed elements
  !0 - ballistic kernel
  !1 - RANDSVD with equal singular values
  !2 - RANDSVD with sigma_k = 1/k
  !3 - RANDSVD with sigma_k = 1/k^2
  !4 - RANDSVD with sigma_k = 1/2^k
  Type(Mtrx) mask !Sparse mask (zeros and ones)
  Type(Mtrx) mat, p1, p2, p3 !Temporary matrices
  Type(Mtrx) u, s, v !SVD factors
  Type(Mtrx) lcurmat, rcurmat !Low rank factors
  Type(Mtrx) curmat !=lcurmat*rcurmat
  Type(SparseRow) KnownSparse !Known elements as a sparse matrix
  Type(Mtrx) matmask !Full version of KnownSparse
  Integer(4) n, r !Size and rank
  Integer(4) curr !Starting rank
  Double precision curtime !Time measurment
  Double precision q !Sparsity
  Double precision dsecnd !Lapack subroutine to measure time
  Character(50) filename !Integer part of the name for output
  Integer(4) svdtime !Number of svds, measuring time
  Integer(4) maxi !Total number of iterations
  Integer(4) i, k
  n = 1000
  print *, 'alpha', alpha
  print *, 'beta', mt
  write(filename,'(i0)') mt
  close(2)
  open(2, file=trim("remake_qrt50_e1.")//trim(filename), position="append")
  !open(2, file=trim("remake_qrt_test.")//trim(filename), position="append")
  do r = 5,25
  
    q = 0.1d0 + 0.015d0*alpha
    svdtime = 50
    
    call p1%random(n,r)
    call p3%random(n,r)
    
    if (mt == 0) then
      call mat%badrandom(n, n, r)
    else if (mt == -1) then
      call mat%init(n, n)
      do i = 1, n
        do k = 1, n
          mat%d(i,k) = (i**(1/3.0d0) + k**(1/3.0d0))**2 * (1.0d0/i + 1.0d0/k)**0.5d0
        end do
      end do
    else if (mt == 1) then
      call p2%init(r,r)
      do i = 1, r
        p2%d(i,i) = 1.0d0
      end do
      mat = p1*p2*(.T.p3)
      call p2%deinit()
    else if (mt == 2) then
      call p2%init(r,r)
      do i = 1, r
        p2%d(i,i) = 1.0d0/i
      end do
      mat = p1*p2*(.T.p3)
      call p2%deinit()
    else if (mt == 3) then
      call p2%init(r,r)
      do i = 1, r
        p2%d(i,i) = 1.0d0/i**2
      end do
      mat = p1*p2*(.T.p3)
      call p2%deinit()
    else if (mt == 4) then
      call p2%init(r,r)
      do i = 1, r
        p2%d(i,i) = 2.0d0/2.0d0**i
      end do
      mat = p1*p2*(.T.p3)
      call p2%deinit()
    end if
    
    curr = 1
    if (mt == 1) then
      curr = r
    end if
    call mask%mask(n, q)
    matmask = mask .dot. mat
    KnownSparse = matmask
    call lcurmat%init(n,curr)
    call rcurmat%init(curr,n)
    
    curtime = dsecnd()
    curtime = dsecnd()
    do i = 1, 1
      call mat%svd(u,s,v)
    end do
    curtime = (dsecnd()-curtime)*svdtime
    print *, 'time_svd', curtime
    
    call reconstruct(KnownSparse, r, lcurmat, rcurmat, curtime, eps, 1, maxi)
    
    curmat = lcurmat * rcurmat
    p1 = mat - curmat
    if (p1%fnorm() > mat%fnorm()) then
      p1 = mat
    end if
    p2 = mask .dot. p1
    print *, -log(p1%fnorm()/mat%fnorm())
    print *, p1%fnorm()/mat%fnorm(), p2%fnorm()/matmask%fnorm()
    write(2, *) q, r, p1%fnorm()/mat%fnorm(), p2%fnorm()/matmask%fnorm(), maxi
    call mask%deinit()
    call mat%deinit()
    call p1%deinit()
    call p2%deinit()
    call matmask%deinit()
    call curmat%deinit()
    call lcurmat%deinit()
    call rcurmat%deinit()
  end do
  close(2)
end
