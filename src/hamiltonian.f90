module hamiltonian

   use data_struct, only: dcsr3

   public:: J2Ham, N0Ham, Ham

contains

   ! create matrix representation for J^2 part
   subroutine J2Ham( A, N, M, c )

      implicit none
      type( dcsr3 ), intent(inout) :: A
      integer, intent(in) :: N, M
      double precision, intent(in) :: c

      real bounds(3)
      integer kmin, kmax, k, nnz, i, nsize

      if( allocated(A%rowind) ) deallocate(A%rowind)
      if( allocated(A%col) ) deallocate(A%col)
      if( allocated(A%val) ) deallocate(A%val)

      ! Even magnetization case
      if( modulo(M, 2) .eq. 0 ) then

         bounds = (/ 0.0, real(M)/2.0, -real(M)/2.0 /)
         kmin   = ceiling( maxval( bounds ) )

         bounds = (/ real(N)/2.0, real(N)-real(M)/2.0, real(N)+real(M)/2.0 /)
         kmax   = floor( minval( bounds ) )

         ! number of nonzero elements
         nsize  = (kmax-kmin+1)
         nnz    = 2*(nsize)-1

         allocate( A%val(nnz), A%col(nnz), A%rowind(nsize+1) )

         A%rowind(1) = 1
         A%m = nsize

         do k=kmin, kmax

            ! update values and columns
            i = A%rowind(k-kmin+1)

            A%val(i) = 1.d0/( 2.d0*dble(N) )*( dble(N) + dble(M)**2 + &
                       ( dble(N) - 2.d0*dble(k) )*( 2.d0*dble(N) + 1.d0 - &
                       2.d0*( dble(N) - 2.d0*dble(k) ) ) )
            A%col(i) = k-kmin+1

            A%val(i) = c * A%val(i)

            if( k<kmax ) then
               A%val(i+1) = 1.d0/dble(N)*SQRT( ( dble(k) + 1.d0 + dble(M)/2.d0 )*&
                            ( dble(k) + 1.d0 - dble(M)/2.d0 )*( dble(N) - 2.d0*dble(k) )*&
                            ( dble(N) - 2.d0*dble(k) - 1.d0 ) )
               A%col(i+1) = k-kmin+2

               A%val(i+1) = c * A%val(i+1)

            endif

            ! update csr rowindex
            if( k< kmax ) then
               i = k-kmin+2
               A%rowind(i) = A%rowind(i-1) + 2
            else
               i = k-kmin+2
               A%rowind(i) = A%rowind(i-1) + 1
            endif

         enddo

      else

         bounds = (/ 0.0, real(M+1)/2.0, -real(M-1)/2.0 /)
         kmin   = ceiling( maxval( bounds ) )

         bounds = (/ real(N-1)/2.0, real(N)-real(M+1)/2.0, real(N)+real(M-1)/2.0 /)
         kmax   = floor( minval( bounds ) )

         ! number of nonzero elements
         nsize  = (kmax-kmin+1)
         nnz    = 2*(nsize)-1

         allocate( A%rowind(nsize+1), A%col(nnz), A%val(nnz) )

         A%rowind(1) = 1
         A%m = nsize

         do k=kmin, kmax
            
            ! update values and columns
            i = A%rowind(k-kmin+1)

            A%val(i) = 1.d0/( 2.d0*dble(N) )*( dble(N) + dble(M)**2 + &
                       ( dble(N) - 2.d0*dble(k) - 1.d0 )*( 2.d0*dble(N) + 1.d0 - &
                       2.d0*( real(N) - 2.d0*dble(k) - 1.d0 ) ) )            
            A%col(i) = k-kmin+1

            A%val(i) = c * A%val(i)

            if( k<kmax ) then
               A%val(i+1) = 1.d0/dble(N)*SQRT( ( dble(k) + 1.5d0 + dble(M)/2.d0 )*&
                            ( dble(k) + 1.5d0 - dble(M)/2.d0 )*( dble(N) - 2.d0*dble(k) - 1.d0 )*&
                            ( dble(N) - 2.d0*dble(k) - 2.d0 ) )
               A%col(i+1) = k-kmin+2

               A%val(i+1) = c * A%val(i+1)

            endif

            ! update csr rowindex
            if( k< kmax ) then
               i = k-kmin+2
               A%rowind(i) = A%rowind(i-1) + 2
            else
               i = k-kmin+2
               A%rowind(i) = A%rowind(i-1) + 1
            endif

         enddo

      endif

   end subroutine

   ! create N0 matrix representation
   subroutine N0Ham( A, N, M, q )

      implicit none
      type( dcsr3 ), intent(inout) :: A
      integer, intent(in) :: N, M
      double precision, intent(in) :: q

      real bounds(3)
      integer kmin, kmax, k, nnz, i, nsize

      if( allocated(A%rowind) ) deallocate(A%rowind)
      if( allocated(A%col) ) deallocate(A%col)
      if( allocated(A%val) ) deallocate(A%val)

      ! Even magnetization case
      if( modulo(M, 2) .eq. 0 ) then

         bounds = (/ 0.0, real(M)/2.0, -real(M)/2.0 /)
         kmin   = ceiling( maxval( bounds ) )

         bounds = (/ real(N)/2.0, real(N)-real(M)/2.0, real(N)+real(M)/2.0 /)
         kmax   = floor( minval( bounds ) )

         ! number of nonzero elements
         nsize  = (kmax-kmin+1)
         nnz    = nsize

         allocate( A%rowind(nsize+1), A%col(nnz), A%val(nnz) )

         A%rowind(1) = 1
         A%m = nsize

         do k=kmin, kmax

            i = k-kmin+1
            A%rowind(i+1) = A%rowind(i) + 1
            A%col(i)      = i
            A%val(i)      = q * 2.d0*dble(k)

         enddo

      else

         bounds = (/ 0.0, real(M+1)/2.0, -real(M-1)/2.0 /)
         kmin   = ceiling( maxval( bounds ) )

         bounds = (/ real(N-1)/2.0, real(N)-real(M+1)/2.0, real(N)+real(M-1)/2.0 /)
         kmax   = floor( minval( bounds ) )

         ! number of nonzero elements
         nsize  = (kmax-kmin+1)
         nnz    = nsize

         allocate( A%rowind(nsize+1), A%col(nnz), A%val(nnz) )

         A%rowind(1) = 1
         A%m = nsize

         do k=kmin, kmax

            i = k-kmin+1
            A%rowind(i+1) = A%rowind(i) + 1
            A%col(i)      = i
            A%val(i)      = q * ( 2.d0*dble(k) + 1.d0 )

         enddo

      endif

   end subroutine

   !> Information about c and q is stored in 
   !> matrixes J2 and N0
   subroutine Ham( A, J2, N0 )

      implicit none
      type( dcsr3 ), intent(inout) :: A
      type( dcsr3 ), intent(in) :: J2, N0

      character trans
      integer request
      integer m, n, nzmax, info, sort, nnz

      trans = 'N'
      m = J2%m
      n = maxval(J2%col)

      if( allocated(A%rowind) ) then
         if( size(A%rowind,1) .ne. m+1 ) then
            deallocate(A%rowind)
            allocate( A%rowind(m+1) )
            A%m = m
         endif
      else
         allocate( A%rowind(m+1) )
         A%m = m
      endif

      ! double call to mkl_dcsradd
      request=1
      call mkl_dcsradd( trans, request, sort, m, n,&
                        J2%val, J2%col, J2%rowind,&
                        1.d0,&
                        N0%val, N0%col, N0%rowind,&
                        A%val, A%col, A%rowind,&
                        nzmax, info )

      nnz = A%rowind(m+1)-1   

      if( allocated(A%col) ) then
         if( size(A%col,1) .ne. m+1 ) then
            deallocate(A%col)
            allocate( A%col(nnz) )
         endif
      else
         allocate( A%col(nnz) )
      endif

      if( allocated(A%val) ) then
         if( size(A%val,1) .ne. m+1 ) then
            deallocate(A%val)
            allocate( A%val(nnz) )
         endif
      else
         allocate( A%val(nnz) )
      endif

      request=2
      call mkl_dcsradd( trans, request, sort, m, n,&
                        J2%val, J2%col, J2%rowind,&
                        1.d0,&
                        N0%val, N0%col, N0%rowind,&
                        A%val, A%col, A%rowind,&
                        nzmax, info )

   end subroutine

end module
