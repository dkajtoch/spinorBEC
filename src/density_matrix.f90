module density_matrix

   use data_struct

   private
   public :: thermal_Mrho, thermal_rho

contains

   function thermal_Mrho( N, M, c, q, beta ) result( res )

      ! Input:
      !
      !  N - number of particles
      !  M - magnetization
      !  c - parameter multiplying J^2 term
      !  q - parameter multiplying (N+ + N-)
      !  H = c*J^2 + q*(N_{+} + N_{-})
      !  beta - temperature 
      !
      ! Output:
      !
      !  res - Mrho data type

      !     res%eigval(i) - exp( -beta * e(i) )/Z
      !     eigenvalues are stored in ascending order
      !     res%eigvec - matrix with eigenvectors (column wise)

      implicit none
      integer, intent(in) :: N, M
      double precision, intent(in) :: c, q, beta
      type( Mrho ) :: res

      ! local variables
      double precision Zstat, E0
      ! for lapack routine
      double precision, allocatable ::  e(:), work(:)
      integer ldz, nn, info
      character*1 compz

      ! for constrution of the Hamiltonian
      integer kmin, kmax, k, i

      ! -----------------------------------------
      ! ----- Construct Hamiltonian matrix ------
      ! -----------------------------------------

      ! Even magnetization case

      if( modulo(M, 2) .eq. 0 ) then

         call EvenBounds( N, M, kmin, kmax )

         ! size of the matrix 
         nn = (kmax - kmin + 1)
         
         ! allocate resources
         if( allocated( res%eigval) ) then
           deallocate( res%eigval )
           allocate( res%eigval(nn) )
         else 
           allocate( res%eigval(nn) )
         endif

         if( allocated( e ) ) then
            deallocate(e)
            allocate( e(nn-1) )
         else
            allocate( e(nn-1) )
         endif

         do k=kmin, kmax

            ! array index
            i = k - kmin + 1

            ! diagonal part
            res%eigval(i) = &

                   c * 1.d0/( 2.d0*dble(N) )*( dble(N) + dble(M)**2 + &
                       ( dble(N) - 2.d0*dble(k) )*( 2.d0*dble(N) + 1.d0 - &
                       2.d0*( dble(N) - 2.d0*dble(k) ) ) ) + &
                   q * (2.d0 * dble(k) - dble(N))
                   !q * (2.d0 * dble(k))

            ! offdiagonal part
            if( k < kmax )  then

               e(i) = c * 1.d0/dble(N)*SQRT( ( dble(k) + 1.d0 + dble(M)/2.d0 )*&
                          ( dble(k) + 1.d0 - dble(M)/2.d0 )*( dble(N) - 2.d0*dble(k) )*&
                          ( dble(N) - 2.d0*dble(k) - 1.d0 ) )

            endif
         
         enddo

      ! Odd magnetization
      else

         call OddBounds( N, M, kmin, kmax )

         ! number of nonzero elements
         nn  = (kmax-kmin+1)

         ! allocate resources
         if( allocated( res%eigval) ) then
           deallocate( res%eigval )
           allocate( res%eigval(nn) )
         else
           allocate( res%eigval(nn) )
         endif

         if( allocated( e ) ) then
            deallocate(e)
            allocate( e(nn-1) )
         else
            allocate( e(nn-1) )
         endif

         do k=kmin, kmax
            
            ! array index
            i = k - kmin + 1

            ! diagonal part
            res%eigval(i) = &
               
                   c * 1.d0/( 2.d0*dble(N) )*( dble(N) + dble(M)**2 + &
                       ( dble(N) - 2.d0*dble(k) - 1.d0 )*( 2.d0*dble(N) + 1.d0 - &
                       2.d0*( real(N) - 2.d0*dble(k) - 1.d0 ) ) ) + &           
                   q * (2.d0*dble(k) + 1.d0 - dble(N))
                   !q * (2.d0*dble(k))

            ! off-diagonal terms
            if( k < kmax ) then

               e(i) = c * 1.d0/dble(N)*SQRT( ( dble(k) + 1.5d0 + dble(M)/2.d0 )*&
                          ( dble(k) + 1.5d0 - dble(M)/2.d0 )*( dble(N) - 2.d0*dble(k) - 1.d0 )*&
                          ( dble(N) - 2.d0*dble(k) - 2.d0 ) )

            endif

         enddo

      endif

      ! ------------------------------------------
      ! --- Diagonalization of the Hamiltonian ---
      ! ------------------------------------------

      compz = 'I'
      ldz   = nn

      ! allocate workspace
      if( allocated( res%eigvec ) ) then
         deallocate( res%eigvec )
         allocate( res%eigvec(nn,nn) )
      else
         allocate( res%eigvec(nn,nn) )
      endif

      if( allocated( work ) ) then
         deallocate( work )
         allocate( work(2*nn-2) )
      else
         allocate( work(2*nn-2) )
      endif

      call DSTEQR( compz, nn, res%eigval, e, res%eigvec, ldz, work, info )

      if( info .ne. 0 ) then
         if( info < 0 ) then
            write(*,'(A,I2.2,A)') 'the ', abs(info), '-th parameter has an illegeal value'
         else
            write(*,'(A,I8.8,A,I8.8,A)')  'the algorithm failed to find all eigenvalues after ' &
               , 30*nn, ' iterations: ', info, ' off-diagonal elements have not converged to zero'
         endif
      endif

      ! return thermal probabilities
      ! avoids overflow exception
      ! floating-point underflow replaced by 0
      E0 = MINVAL( res%eigval )
      Zstat = 0.0d0
      do i = 1, nn
         Zstat = Zstat + EXP( -beta * ( res%eigval(i) - E0 ) )
      enddo

      do i = 1, nn
         res%eigval(i) = EXP( -beta * ( res%eigval(i) - E0 ) )/Zstat
      enddo

   end function

   function thermal_rho( N, mu, sigma, c, q, beta ) result(res)

   ! INPUT:
   !
   !     N    - number of particles
   !     Mbar - center of magnetization 
   !     sigma - dispersion of magnetization
   !     c - parameter of the Hamiltonian J^2
   !     q - parameter of the Hamiltonian N_0
   !     beta = 1/(k_B T)
   !
   ! OUTPUT:
   !
   !     array of data type( Mrho )

      implicit none
      integer, intent(in) :: N
      real( kind = dp ) :: mu, sigma
      double precision, intent(in) :: c, q, beta
      type( rho ) :: res

      ! local variables
      real( kind = dp ) :: norm
      integer M
      ! for lapack routine
      double precision, allocatable ::  e(:), work(:)
      integer ldz, nn, info
      character*1 compz

      ! for constrution of the Hamiltonian
      integer kmin, kmax, k, i, midx 

      if( allocated(res%array) ) then
         deallocate( res%array )
         allocate( res%array(0:N) )
      else
         allocate( res%array(0:N) )
      endif

      ! ------------------------------------------------
      !            Loop over magnetizations
      ! ------------------------------------------------
      ! parallel with OpenMP

      ! $OMP PARALLEL DO
      do M = 0, N
         res%array(M) = thermal_Mrho( N, M, c, q, beta )
      enddo
      ! $OMP END PARALLEL DO

      ! -------------------------------------
      ! Gaussian weights for magentization
      ! subspaces M = -N, ..., 0. ,,, N
      ! -------------------------------------
      if( allocated( res%weights) ) then
         deallocate( res%weights )
         allocate( res%weights(-N:N) )
      else
         allocate( res%weights(-N:N) )
      endif

      ! normalized weights
      norm = 0.0_dp
      do M = -N, N
         res%weights(M) = EXP( -0.5_dp * ( (real(M,dp)- mu)/sigma )**2 )
         norm = norm + res%weights(M)
      enddo

      do M = -N, N
         res%weights(M) = res%weights(M)/norm
      enddo

   end function

end module
