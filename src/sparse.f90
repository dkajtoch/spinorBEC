module sparse

   ! use mkl include files
   !use mkllib

   use data_struct, only: dcsr3

   private
   public:: eigen
   
contains

!  | -------------------------------------------------- |
!  | Perform matrix-vector product using:               |
!  | upper-symmetric-sparse-csr matrix and dense vector |
!  | -------------------------------------------------- |
   subroutine dcsr3symv( A, x, y )

      implicit none
      type( dcsr3 ), intent(in) :: A
      double precision, intent(in) :: x(:)
      double precision, intent(out) :: y(:)
         
      character :: uplo
      integer :: m

      ! take number of rows
      m = A%m
      ! choose upper diagonal 
      uplo = 'U'

      call mkl_dcsrsymv( uplo, m, A%val, A%rowind, A%col, &
                         x, y )

   end subroutine dcsr3symv

!  | --------------------------------------------- |
!  | Find eigenvalues and eigenvectors of real     |
!  | symmetric matrix using ARPACK routines        |
!  | --------------------------------------------- |
   subroutine eigen( A, d, z, nev, which, op_ncv, op_tol, op_maxiter)

     implicit none
     type( dcsr3 ), intent(in) :: A
     double precision, intent(inout), allocatable :: d(:)   ! eigenvalues
     double precision, intent(inout), allocatable :: z(:,:) ! eigenvectors (stored in columns)
     integer, intent(in) :: nev ! number of eigenvalues
     character, intent(in) :: which(2) ! 'LA','SA','LM','SM','BE'
     integer, optional, intent(in) :: op_ncv
     integer, optional, intent(in) :: op_maxiter
     double precision, optional, intent(in) :: op_tol

     ! arpack parameters
     character bmat
     integer sigma
     integer ido
     integer info
     integer iparam(11)
     integer ipntr(11)
     integer lworkl
     integer ldv
     integer n
     integer ncv
     double precision tol
     logical, allocatable :: select(:)
     double precision, allocatable :: resid(:)
     double precision, allocatable :: workl(:)
     double precision, allocatable, target :: workd(:)
     double precision, allocatable :: v(:,:)
     double precision, pointer :: x(:) => null()
     double precision, pointer :: y(:) => null()

     ido    = 0           ! reverse communication interface init
     info   = 0           ! info output
     bmat   = 'I'         ! standard eigenvalue problem
     n      = A%m         ! dimension of the eigenpoblem
     ncv    = 2*nev + 1   ! number of Lanczos vectors
     tol    = EPSILON(0.0)! tolerance for convergence

     if( ncv > n ) ncv=n
     if( present(op_ncv) ) ncv=op_ncv
     if( present(op_tol) ) tol=op_tol

     lworkl = ncv*(ncv+8) ! workspace size

     allocate( resid(n) ) 
     allocate( workd(3*n) )
     allocate( workl(lworkl) )
     allocate( v(n,ncv) )
     allocate( select(ncv) )

     ldv    = n           ! leading dimension of v

     ! | -------------- |
     ! | setting iparam |
     ! | -------------- |
     iparam(1) = 1       ! exact shifts strategy
     iparam(3) = 10000   ! maximum number of iterations
     iparam(7) = 1       ! mode - standard eigenvalue

     if( present(op_maxiter) ) iparam(3)=op_maxiter

     ! | -------------------------- |
     ! | reverse communication loop |
     ! | -------------------------- |
     do

        call dsaupd( ido, bmat, n, which, nev, tol, resid, ncv, v, &
                     ldv, iparam, ipntr, workd, workl, lworkl, info )

        if( ido .ne. 1 .and. ido .ne. -1 ) then
           exit
        endif

        ! perform matrix-vector multiplication
        x => workd( ipntr(1):ipntr(1)+n )
        y => workd( ipntr(2):ipntr(2)+n )

        call mkl_dcsrsymv( 'U', n, A%val, A%rowind, A%col, &
                           x, y )

        !call dcsr3symv( A, x, y )

      enddo

      if( associated(x) ) nullify(x)
      if( associated(y) ) nullify(y)

      ! | -------------------- |
      ! | Check for errors     |
      ! | -------------------- |
      if ( info < 0 ) then

         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'EIGEN - Fatal error!'
         write ( *, '(a,i6)' ) 'Error with EIGEN, INFO = ', info
         write ( *, '(a)' ) 'Check documentation in DSAUPD.'

      elseif( info == 1 ) then

         write(*, '(a)' ) ' '
         write(*, '(a)' ) 'ARPACK Warning!'
         write(*, '(a)' ) 'Maximum number of interations taken'
         write(*, '(a,I0)' ) 'Number of converged Ritz values: ', iparam(5)
         write(*, '(a)' ) 'Try changing maximum number of iterations (op_maxiter) or number of Arnoldi vectors (op_ncv)'

      endif 
         

      !
      !  No fatal errors occurred.
      !  Post-Process using DSEUPD.
      !
      !  Computed eigenvalues may be extracted.
      !
      !  Eigenvectors may be also computed now if
      !  desired.  (indicated by rvec = .true.)
      !
      !  The routine DSEUPD now called to do this
      !  post processing
      !

         if( allocated(d) ) then
            if( size(d,1) .ne. nev ) then
               deallocate( d )
               allocate( d(nev) )
            endif
         else
            allocate( d(nev) )
         endif

         if( allocated(z) ) then
            if( size(z,1) .ne. n .or. size(z,2) .ne. nev ) then
               deallocate( z )
               allocate( z(n,nev) )
            endif
         else
            allocate( z(n,nev) )
         endif

         call dseupd ( .true., 'A', select, d, z, ldv, sigma, &
                       bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                       iparam, ipntr, workd, workl, lworkl, info )

         if ( info .ne. 0 ) then

            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'EIGEN - Fatal error!'
            write ( *, '(a,i6)' ) 'Error with DSEUPD, INFO = ', info
            write ( *, '(a)' ) 'Check the documentation of DSEUPD.'

         endif

  end subroutine eigen

end module sparse
