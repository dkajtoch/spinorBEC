module measures

   use sparse, only: eigen
   use data_struct
   use density_matrix, only: thermal_Mrho

   private
   public:: EnergyGap, FisherQ, Mean0, Var0, Fluct, DensityCorrelation, InverseMean0

   interface Mean0
      module procedure :: vec_Mean0
      module procedure :: rhoM_Mean0
      module procedure :: rho_Mean0
   end interface

   interface InverseMean0
      module procedure :: rhoM_InverseMean0
   end interface

   interface FisherQ
      module procedure :: vec_FisherQ
      module procedure :: Mrho_FisherQ
      module procedure :: rho_FisherQ
   end interface

   interface Fluct
      module procedure :: vec_Fluct
      module procedure :: Mrho_Fluct
      module procedure :: rho_Fluct
   end interface

   contains

   function EnergyGap( A, maxiter, ncv )

      implicit none
      type( dcsr3 ), intent(in) :: A
      integer, optional :: maxiter, ncv
      double precision EnergyGap

      ! local variables
      double precision, allocatable :: eigval(:), eigvec(:,:)

      if( present(maxiter) ) then
         if( present(ncv) ) then
            call eigen( A, eigval, eigvec, 2, 'SA', op_maxiter=maxiter, op_ncv=ncv )
         else
            call eigen( A, eigval, eigvec, 2, 'SA', op_maxiter=maxiter )
         endif
      else
         if( present(ncv) ) then
            call eigen( A, eigval, eigvec, 2, 'SA', op_ncv=ncv )
         else
            call eigen( A, eigval, eigvec, 2, 'SA' )
         endif
      endif

      EnergyGap = abs( eigval(2) - eigval(1) )

   end function

   function vec_Mean0( vec, N, M ) result(Mean0)

      implicit none
      double precision, intent(in) :: vec(:)
      integer, intent(in) :: N, M
      double precision :: Mean0

      ! Local variables
      integer kmin, kmax, k, i

      Mean0 = 0.d0

      ! Even magnetization case
      if( modulo(M, 2) .eq. 0 ) then

         call EvenBounds( N, M, kmin, kmax )

         do k=kmin, kmax

            ! index shift
            i = k-kmin+1
            Mean0 = Mean0 + ( dble(N) - 2.d0*dble(k) ) * vec(i)**2

         enddo

      else

         call OddBounds( N, M, kmin, kmax )

         do k=kmin, kmax

            i = k-kmin+1
            Mean0 = Mean0 + ( dble(N) - 2.d0*dble(k) - 1.d0 ) * vec(i)**2

         enddo

      endif

   end function

   function rhoM_Mean0( rhoin, N, M ) result( Mean0 )

      implicit none
      type( Mrho ), intent(in) :: rhoin
      integer, intent(in) :: N, M
      double precision :: Mean0

      ! local variables
      integer idx, i, k ,kmin, kmax

      Mean0 = 0.d0

      ! Even magnetization case
      if( modulo(M,2) == 0 ) then

         call EvenBounds( N, M, kmin, kmax )

         ! loop over eigenvectors
         do idx = 1, size( rhoin%eigval )
            do k = kmin, kmax
               i = k-kmin+1
               Mean0 = Mean0 + rhoin%eigval(idx) * ( dble(N) - 2.d0*dble(k) ) * rhoin%eigvec(i,idx)**2
            enddo
         enddo

      ! Odd magnetization case
      else

         call OddBounds( N, M, kmin, kmax )

         ! loop over eigenvectors
         do idx = 1, size( rhoin%eigval )
            do k = kmin, kmax
               i = k-kmin+1
               Mean0 = Mean0 + rhoin%eigval(idx) * ( dble(N) - 2.d0*dble(k) - 1.d0 ) * rhoin%eigvec(i,idx)**2
            enddo
         enddo

      endif

   end function

   function rho_Mean0( rhoin, N ) result( Mean0 )

      implicit none
      type( rho ), intent(in) :: rhoin
      integer, intent(in) :: N
      double precision :: Mean0

      ! local variables
      integer M, i, idx, k, kmin, kmax

      Mean0 = 0.d0

      ! loop over magnetizations
      do M = -N, N

         ! even magnetization case
         if( modulo(M,2) == 0 ) then

            call EvenBounds( N, M, kmin, kmax )

            ! loop over eigenvectors
            do idx = 1, size( rhoin%array(abs(M))%eigval )
               do k = kmin, kmax
                  i = k-kmin+1
                  Mean0 = Mean0 + rhoin%array(abs(M))%eigval(idx) * ( dble(N) - 2.d0 * dble(k) ) * &
                                  rhoin%array(abs(M))%eigvec(i,idx)**2
               enddo
            enddo

         ! Odd magnetization case
         else

            call OddBounds( N, M, kmin, kmax )

            ! loop over eigenvectors
            do idx = 1, size( rhoin%array(abs(M))%eigval )
               do k = kmin, kmax
                  i = k-kmin+1
                  Mean0 = Mean0 + rhoin%array(abs(M))%eigval(idx) * ( dble(N) - 2.d0 * dble(k) - 1.d0 ) * &
                                  rhoin%array(abs(M))%eigvec(i,idx)**2
               enddo
            enddo

         endif

      enddo

   end function

   subroutine rhoM_InverseMean0( rhoin, N, M, c, beta, n0, q )

      implicit none
      type( Mrho ), intent(inout) :: rhoin
      integer, intent(in) :: N, M
      double precision, intent(in) :: c, beta, n0
      double precision, intent(out) :: q

      ! local variables
      double precision :: a, b, d, fa, fb, fd
      double precision :: n0max, n0min
      integer :: kmin, kmax
      double precision, parameter :: eps = 1.d-08

      ! maximal and minimal values
      if( modulo(M,2) == 0 ) then
         call EvenBounds( N, M, kmin, kmax )
         n0max = 1.d0 - 2.d0 * dble(kmin)/dble(N)
         n0min = 1.d0 - 2.d0 * dble(kmax)/dble(N)
      else
         call OddBounds( N, M, kmin, kmax )
         n0max = 1.d0 - 2.d0 * dble(kmin)/dble(N) - 1.d0/dble(N)
         n0min = 1.d0 - 2.d0 * dble(kmax)/dble(N) - 1.d0/dble(N)
      endif

      if( n0 > n0max .or. n0 < n0min ) then
         print *, 'Incorrect n0 value'
         q = 0.d0
      else
         ! search for initial conditions
         ! fa * fb < 0.d0
         a = 1.d+03
         b = -1.d+03

         rhoin = thermal_Mrho( N, M, c, a, beta )
         fa    = Mean0( rhoin, N, M ) - N * n0
         rhoin = thermal_Mrho( N, M, c, b, beta )
         fb    = Mean0( rhoin, N, M ) - N * n0

         do while( fa * fb > 0.d0 )
            ! bad starting conditions
            a = 2.d0 * a
            b = 2.d0 * b

            rhoin = thermal_Mrho( N, M, c, a, beta )
            fa    = Mean0( rhoin, N, M ) - N * n0
            rhoin = thermal_Mrho( N, M, c, b, beta )
            fb    = Mean0( rhoin, N, M ) - N * n0
         enddo

         if( fa == 0.d0 ) then
            q = a
         elseif( fb == 0.d0 ) then
            q = b
         else
            ! search for 0 with bisection
            do while( .true. )
               d = (a+b)/2.d0
               rhoin = thermal_Mrho( N, M, c, d, beta )
               fd    = Mean0( rhoin, N, M ) - n0 * N

               if( fa * fd < 0.d0 ) then
                  b  = d
                  fb = fd
               else
                  a  = d
                  fa = fd
               endif

               ! stop condition
               if( abs(fd) < eps) then
                  q = d
                  exit
               endif
            enddo
         endif
      endif

   end subroutine

   function Var0( vec, N, M )

      implicit none
      double precision, intent(in) :: vec(:)
      integer, intent(in) :: N, M
      double precision :: Var0

      ! Local variables
      double precision sumn, sumn2
      integer kmin, kmax, k, i

      sumn = 0.d0; sumn2 = 0.d0

      ! Even magnetization case
      if( modulo(M, 2) .eq. 0 ) then

         call EvenBounds( N, M, kmin, kmax )

         do k=kmin, kmax

            ! index shift
            i = k-kmin+1
            sumn  = sumn  + ( dble(N) - 2.d0*dble(k) ) * vec(i)**2
            sumn2 = sumn2 + ( dble(N) - 2.d0*dble(k) )**2 * vec(i)**2

         enddo

      else

         call OddBounds( N, M, kmin, kmax )

         do k=kmin, kmax

            i = k-kmin+1
            sumn  = sumn  + ( dble(N) - 2.d0*dble(k) - 1.d0 ) * vec(i)**2
            sumn2 = sumn2 + ( dble(N) - 2.d0*dble(k) - 1.d0 )**2 * vec(i)**2

         enddo

      endif

      Var0 = sumn2 - sumn**2

   end function

   function DensityCorrelation( vec, N, M ) result( vout )

      implicit none
      double precision, intent(in) :: vec(:)
      integer, intent(in) :: N, M
      double precision :: vout(3,3)

      ! Local variables
      double precision sumP, sum0, sumM, sumPP, sumP0, sumPM, sum00, sum0M, sumMM
      integer kmin, kmax, k, i

      sumP = 0.d0
      sum0 = 0.d0
      sumM = 0.d0
      sumPP = 0.d0
      sumP0 = 0.d0
      sumPM = 0.d0
      sum00 = 0.d0
      sum0M = 0.d0
      sumMM = 0.d0

      ! Even magnetization case
      if( modulo(M, 2) .eq. 0 ) then

         call EvenBounds( N, M, kmin, kmax )

         do k=kmin, kmax

            ! index shift
            i = k-kmin+1
            sumP = sumP + ( dble(k) + dble(M/2) ) * vec(i)**2
            sum0 = sum0 + ( dble(N) - 2.d0*dble(k) ) * vec(i)**2
            sumM = sumM + ( dble(k) - dble(M/2) ) * vec(i)**2
            sumPP = sumPP + ( dble(k) + dble(M/2) )**2 * vec(i)**2
            sumP0 = sumP0 + ( dble(k) + dble(M/2) ) * ( dble(N) - 2.d0*dble(k) ) * vec(i)**2
            sumPM = sumPM + ( dble(k) + dble(M/2) ) * ( dble(k) - dble(M/2) ) * vec(i)**2
            sum00 = sum00 + ( dble(N) - 2.d0*dble(k) )**2 * vec(i)**2
            sum0M = sum0M + ( dble(N) - 2.d0*dble(k) ) * ( dble(k) - dble(M/2) ) * vec(i)**2
            sumMM = sumMM + ( dble(k) - dble(M/2) )**2 * vec(i)**2

         enddo

      else

         call OddBounds( N, M, kmin, kmax )

         do k=kmin, kmax

            i = k-kmin+1
            sumP = sumP + ( dble(k) + dble((M+1)/2) ) * vec(i)**2
            sum0 = sum0 + ( dble(N) - 2.d0*dble(k) - 1.d0 ) * vec(i)**2
            sumM = sumM + ( dble(k) - dble((M-1)/2) ) * vec(i)**2
            sumPP = sumPP + ( dble(k) + dble((M+1)/2) )**2 * vec(i)**2
            sumP0 = sumP0 + ( dble(k) + dble((M+1)/2) ) * ( dble(N) - 2.d0*dble(k) - 1.d0 ) * vec(i)**2
            sumPM = sumPM + ( dble(k) + dble((M+1)/2) ) * ( dble(k) - dble((M-1)/2) ) * vec(i)**2
            sum00 = sum00 + ( dble(N) - 2.d0*dble(k) - 1.d0 )**2 * vec(i)**2
            sum0M = sum0M + ( dble(N) - 2.d0*dble(k) - 1.d0 ) * ( dble(k) - dble((M-1)/2) ) * vec(i)**2
            sumMM = sumMM + ( dble(k) - dble((M-1)/2) )**2 * vec(i)**2

         enddo

      endif

      vout(1,1) = sumPP - sumP**2
      vout(1,2) = sumP0 - sumP*sum0; vout(2,1) = sumP0 - sumP*sum0
      vout(1,3) = sumPM - sumP*sumM; vout(3,1) = sumPM - sumP*sumM
      vout(2,2) = sum00 - sum0**2
      vout(2,3) = sum0M - sum0*sumM; vout(3,2) = sum0M - sum0*sumM
      vout(3,3) = sumMM - sumM**2

   end function

   ! ----------------------------------------------------
   ! Return covariance matrix for pure states defined 
   ! within subspace of fixed magnetization
   ! ----------------------------------------------------
   function vec_Covariance( vec, N, M ) result (cov)

      implicit none
      double precision, intent(in) :: vec(:) 
      integer, intent(in) :: N, M
      double precision :: cov(8,8)

      ! local variables
      double precision sumP, sum0, sumM, sumPM, sumP0, sum0M, sumPM00
      double precision sumPP, sum00, sumMM
      integer k, kmin, kmax, i

      sumP  = 0.d0
      sum0  = 0.d0
      sumM  = 0.d0
      sumPP = 0.d0
      sumMM = 0.d0
      sum00 = 0.d0
      sumPM = 0.d0
      sumP0 = 0.d0
      sum0M = 0.d0
      sumPM00 = 0.d0

      ! Even magnetization case
      if( modulo(M, 2) .eq. 0 ) then

         call EvenBounds( N, M, kmin, kmax )

         do k=kmin, kmax

            ! index shift
            i = k-kmin+1

            sumP  = sumP + ( dble(k + M/2) ) * vec(i)**2
            sum0  = sum0 + ( dble(N - 2*k) ) * vec(i)**2
            sumM  = sumM + ( dble(k - M/2) ) * vec(i)**2
            sumPP = sumPP + ( dble(k + M/2)**2 ) * vec(i)**2
            sum00 = sum00 + ( dble(N - 2*k)**2 ) * vec(i)**2
            sumMM = sumMM + ( dble(k - M/2)**2 ) * vec(i)**2
            sumPM = sumPM + ( dble(k + M/2) * dble(k - M/2) ) * vec(i)**2
            sumP0 = sumP0 + ( dble(k + M/2) * dble(N - 2*k) ) * vec(i)**2
            sum0M = sum0M + ( dble(N - 2*k) * dble(k - M/2) ) * vec(i)**2

            if( k < kmax )  then
               sumPM00 = sumPM00 + SQRT( ( dble(k+1+M/2) ) * ( dble(k+1-M/2) ) *&
                         ( dble(N-2*k) ) * ( dble(N-2*k-1) ) ) * vec(i+1) * vec(i)
            endif

         enddo

      else

         call OddBounds( N, M, kmin, kmax )

         do k=kmin, kmax

            i = k-kmin+1

            sumP  = sumP + ( dble(k + (M-1)/2) ) * vec(i)**2
            sum0  = sum0 + ( dble(N - 2*k - 1) ) * vec(i)**2
            sumM  = sumM + ( dble(k - (M+1)/2) ) * vec(i)**2
            sumPP = sumPP + ( dble(k + (M-1)/2)**2 ) * vec(i)**2
            sum00 = sum00 + ( dble(N - 2*k - 1)**2 ) * vec(i)**2
            sumMM = sumMM + ( dble(k - (M+1)/2)**2 ) * vec(i)**2
            sumPM = sumPM + ( dble(k + (M-1)/2) * dble(k - (M+1)/2) ) * vec(i)**2
            sumP0 = sumP0 + ( dble(k + (M-1)/2) * dble(N - 2*k - 1) ) * vec(i)**2
            sum0M = sum0M + ( dble(N - 2*k - 1) * dble(k - (M+1)/2) ) * vec(i)**2

            if( k < kmax )  then
               sumPM00 = sumPM00 + SQRT( ( dble(k+1+(M-1)/2) ) * ( dble(k+1-(M+1)/2) ) *&
                         ( dble(N-2*k-1) ) * ( dble(N-2*k-2) ) ) * vec(i+1) * vec(i)
            endif

         enddo

      endif

      ! put values into the covariance matrix
      cov = 0.d0

      ! Lam = {jx, jy, jz, qxy, qyz, qzx, dxy, y}
      ! ---------------------------------------------------
      ! (1,1)
      cov(1,1) = 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M + 4.d0*sumPM00 + 2.d0*sumP0 )
      ! (1,2)
      ! (1,3)
      ! (1,4)
      ! (1,5)
      ! (1,6)
      cov(1,6) = 0.5d0 * ( -sumM + sumP - 2.d0*sum0M + 2.d0*sumP0 )
      ! (1,7)
      ! (1,8)
      ! ---------------------------------------------------
      ! (2,2)
      cov(2,2) = 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M + 4.d0*sumPM00 + 2.d0*sumP0 )
      ! (2,3)
      ! (2,4)
      ! (2,5)
      cov(2,5) = 0.5d0 * ( -sumM + sumP - 2.d0*sum0M + 2.d0*sumP0 )
      ! (2,6)
      ! (2,7)
      ! (2,8)
      ! ---------------------------------------------------
      ! (3,3)
      cov(3,3) = sumPP + sumMM - 2.d0*sumPM - (sumP - sumM)**2
      ! (3,4)
      ! (3,5)
      ! (3,6)
      ! (3,7)
      ! (3,8)
      cov(3,8) = 1.d0/sqrt(3.d0)*( sumPP - sumMM - 2.d0*sumP0 + 2.d0*sum0M )-&
               1.d0/sqrt(3.d0)*(sumP - sumM)*(sumP - 2.d0*sum0 + sumM)
      ! ---------------------------------------------------
      ! (4,4)
      cov(4,4) = sumM + sumP + 2.d0*sumPM
      ! (4,5)
      ! (4,6)
      ! (4,7)
      ! (4,8)
      ! ---------------------------------------------------
      ! (5,5)
      cov(5,5) = 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M - 4.d0*sumPM00 + 2.d0*sumP0 )
      ! (5,6)
      ! (5,7)
      ! (5,8)
      ! ---------------------------------------------------
      ! (6,6)
      cov(6,6) = 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M - 4.d0*sumPM00 + 2.d0*sumP0 )
      ! (6,7)
      ! (6,8)
      ! ---------------------------------------------------
      ! (7,7)
      cov(7,7) = sumM  + sumP + 2.d0*sumPM
      ! (7,8)
      ! ---------------------------------------------------
      ! (8,8)
      cov(8,8) = 1.d0/3.d0*( sumPP + 4.d0*sum00 + sumMM - 4.d0*sumP0 - 4.d0*sum0M + 2.d0*sumPM )-&
               1.d0/3.d0*(sumP - 2.d0*sum0 + sumM)**2

   end function

   ! ---------------------------------------------------
   ! Optimal quantum Fisher information for pure state
   ! defined within fixed magnetization subspace
   ! --------------------------------------------------
   function vec_FisherQ( vec, N, M ) result(res)

      implicit none
      double precision, intent(in) :: vec(:) 
      integer, intent(in) :: N, M
      type( fis ) :: res

      ! for lapack diagonalization
      integer, parameter :: nn = 8
      integer, parameter :: lda = 8
      integer info, lwork
      double precision a(lda, nn), w(nn)
      double precision :: work(100)

      a = vec_Covariance( vec, N, M )

      ! query the optimal workspace
      lwork = -1
      call DSYEV( 'V', 'U', nn, a, lda, w, work, lwork, info )
      lwork = int( work(1) )

      ! solve eigenvalue
      call DSYEV( 'V', 'U', nn, a, lda, w, work, lwork, info )

      res%val = w(8)/dble(N)**2
      res%dir = a(:,8)

   end function

   ! -----------------------------------------------
   ! Covariance matrix for mixed state defined 
   ! within space of fixed magentization
   ! -----------------------------------------------
   function Mrho_Covariance( rhoin, N, M ) result(cov)

      implicit none
      type( Mrho ), intent(in) :: rhoin 
      integer, intent(in) :: N, M
      real( kind=dp ) :: cov(8,8)

      ! local variables
      double precision sumP, sum0, sumM, sumPM, sumP0, sum0M, sumPM00
      double precision sumPP, sum00, sumMM
      double precision sumJz2, sumY2, sumJzY, sumJz, sumY
      real( kind=dp ) prob, probt
      integer k, kmin, kmax, i, idx, idxi, idxj

      ! flush covariance matrix
      cov(:,:) = 0.0_dp

      ! ====================================================================
      !                     Even magnetization case 
      ! ====================================================================
      if( modulo(M, 2) .eq. 0 ) then

         call EvenBounds( N, M, kmin, kmax )

         ! ===================================
         ! Diagonal part of covariance matrix
         ! ===================================
         do idx = 1, size( rhoin%eigval )

            prob = rhoin%eigval(idx)

            if( prob > 0.0_dp ) then

               sumP  = 0.d0
               sum0  = 0.d0
               sumM  = 0.d0
               sumPP = 0.d0
               sumMM = 0.d0
               sum00 = 0.d0
               sumPM = 0.d0
               sumP0 = 0.d0
               sum0M = 0.d0
               sumPM00 = 0.d0


               do k=kmin, kmax

                  ! index shift
                  i = k-kmin+1

                  sumP  = sumP + ( dble(k + M/2) ) * rhoin%eigvec(i,idx)**2
                  sum0  = sum0 + ( dble(N - 2*k) ) * rhoin%eigvec(i,idx)**2
                  sumM  = sumM + ( dble(k - M/2) ) * rhoin%eigvec(i,idx)**2
                  sumPP = sumPP + ( dble(k + M/2)**2 ) * rhoin%eigvec(i,idx)**2
                  sum00 = sum00 + ( dble(N - 2*k)**2 ) * rhoin%eigvec(i,idx)**2
                  sumMM = sumMM + ( dble(k - M/2)**2 ) * rhoin%eigvec(i,idx)**2
                  sumPM = sumPM + ( dble(k + M/2) * dble(k - M/2) ) * rhoin%eigvec(i,idx)**2
                  sumP0 = sumP0 + ( dble(k + M/2) * dble(N - 2*k) ) * rhoin%eigvec(i,idx)**2
                  sum0M = sum0M + ( dble(N - 2*k) * dble(k - M/2) ) * rhoin%eigvec(i,idx)**2

                  if( k < kmax )  then
                     sumPM00 = sumPM00 + SQRT( ( dble(k+1+M/2) ) * ( dble(k+1-M/2) ) *&
                               ( dble(N-2*k) ) * ( dble(N-2*k-1) ) ) * rhoin%eigvec(i+1,idx) * rhoin%eigvec(i,idx)
                  endif

               enddo

               ! put values into the covariance matrix
               ! Lam = {jx, jy, jz, qxy, qyz, qzx, dxy, y}
               ! ---------------------------------------------------
               cov(1,1) = cov(1,1) + prob * 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M + 4.d0*sumPM00 + 2.d0*sumP0 )
               cov(1,6) = cov(1,6) + prob * 0.5d0 * ( -sumM + sumP - 2.d0*sum0M + 2.d0*sumP0 )
               cov(2,2) = cov(2,2) + prob * 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M + 4.d0*sumPM00 + 2.d0*sumP0 )
               cov(2,5) = cov(2,5) + prob * 0.5d0 * ( -sumM + sumP - 2.d0*sum0M + 2.d0*sumP0 )
               cov(3,3) = cov(3,3) + prob * ( sumPP + sumMM - 2.d0*sumPM - (sumP - sumM)**2 )
               cov(3,8) = cov(3,8) + prob * ( 1.d0/sqrt(3.d0)*( sumPP - sumMM - 2.d0*sumP0 + 2.d0*sum0M )-&
                        1.d0/sqrt(3.d0)*(sumP - sumM)*(sumP - 2.d0*sum0 + sumM) )
               cov(4,4) = cov(4,4) + prob * ( sumM + sumP + 2.d0*sumPM )
               cov(5,5) = cov(5,5) + prob * 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M - 4.d0*sumPM00 + 2.d0*sumP0 )
               cov(6,6) = cov(6,6) + prob * 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M - 4.d0*sumPM00 + 2.d0*sumP0 )
               cov(7,7) = cov(7,7) + prob * ( sumM  + sumP + 2.d0*sumPM )
               cov(8,8) = cov(8,8) + prob * ( 1.d0/3.d0*( sumPP + 4.d0*sum00 + sumMM - 4.d0*sumP0 - 4.d0*sum0M + 2.d0*sumPM )-&
                        1.d0/3.d0*(sumP - 2.d0*sum0 + sumM)**2 )

            endif

         enddo

         ! =======================================
         ! Off-diagonal part of covariance matrix
         ! =======================================
         sumJz2 = 0.0d0
         sumY2  = 0.0d0
         sumJzY = 0.0d0

         do idxi = 1, size( rhoin%eigval )

            do idxj = idxi+1, size( rhoin%eigval )

               prob  = rhoin%eigval(idxi) + rhoin%eigval(idxj)

               if( prob > 0.d0 ) then

                  probt = rhoin%eigval(idxi) * rhoin%eigval(idxj) / prob

                  sumJz = 0.d0
                  sumY  = 0.d0

                  do k=kmin, kmax

                     ! index shift
                     i = k-kmin+1

                     sumJz = sumJz + dble(M) * rhoin%eigvec(i,idxi) * rhoin%eigvec(i,idxj)
                     sumY  = sumY + 2.d0/sqrt(3.d0) * dble(3*k-N) * rhoin%eigvec(i,idxi) * rhoin%eigvec(i,idxj)

                  enddo

                  sumJz2 = sumJz2 + 4.d0 * probt * ( sumJz * sumJz )
                  sumY2  = sumY2  + 4.d0 * probt * ( sumY * sumY )
                  sumJzY = sumJzY + 4.d0 * probt * ( sumJz * sumY )

               endif

            enddo
         enddo

         cov(3,3) = cov(3,3) - sumJz2
         cov(8,8) = cov(8,8) - sumY2
         cov(3,8) = cov(3,8) - sumJzY

      ! =======================================================================
      !                      Odd Magnetization case
      ! =======================================================================
      else

         call OddBounds( N, M, kmin, kmax )

         ! ===================================
         ! Diagonal part of covariance matrix
         ! ===================================
         do idx = 1, size( rhoin%eigval )

            prob = rhoin%eigval(idx)

            if( prob > 0.0_dp ) then

               sumP  = 0.d0
               sum0  = 0.d0
               sumM  = 0.d0
               sumPP = 0.d0
               sumMM = 0.d0
               sum00 = 0.d0
               sumPM = 0.d0
               sumP0 = 0.d0
               sum0M = 0.d0
               sumPM00 = 0.d0

               do k=kmin, kmax

                  i = k-kmin+1

                  sumP  = sumP + ( dble(k + (M-1)/2) ) * rhoin%eigvec(i,idx)**2
                  sum0  = sum0 + ( dble(N - 2*k - 1) ) * rhoin%eigvec(i,idx)**2
                  sumM  = sumM + ( dble(k - (M+1)/2) ) * rhoin%eigvec(i,idx)**2
                  sumPP = sumPP + ( dble(k + (M-1)/2)**2 ) * rhoin%eigvec(i,idx)**2
                  sum00 = sum00 + ( dble(N - 2*k - 1)**2 ) * rhoin%eigvec(i,idx)**2
                  sumMM = sumMM + ( dble(k - (M+1)/2)**2 ) * rhoin%eigvec(i,idx)**2
                  sumPM = sumPM + ( dble(k + (M-1)/2) * dble(k - (M+1)/2) ) * rhoin%eigvec(i,idx)**2
                  sumP0 = sumP0 + ( dble(k + (M-1)/2) * dble(N - 2*k - 1) ) * rhoin%eigvec(i,idx)**2
                  sum0M = sum0M + ( dble(N - 2*k - 1) * dble(k - (M+1)/2) ) * rhoin%eigvec(i,idx)**2

                  if( k < kmax )  then
                     sumPM00 = sumPM00 + SQRT( ( dble(k+1+(M-1)/2) ) * ( dble(k+1-(M+1)/2) ) *&
                               ( dble(N-2*k-1) ) * ( dble(N-2*k-2) ) ) * rhoin%eigvec(i+1,idx) * rhoin%eigvec(i,idx)
                  endif

               enddo

               ! put values into the covariance matrix
               ! Lam = {jx, jy, jz, qxy, qyz, qzx, dxy, y}
               ! ---------------------------------------------------
               cov(1,1) = cov(1,1) + prob * 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M + 4.d0*sumPM00 + 2.d0*sumP0 )
               cov(1,6) = cov(1,6) + prob * 0.5d0 * ( -sumM + sumP - 2.d0*sum0M + 2.d0*sumP0 )
               cov(2,2) = cov(2,2) + prob * 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M + 4.d0*sumPM00 + 2.d0*sumP0 )
               cov(2,5) = cov(2,5) + prob * 0.5d0 * ( -sumM + sumP - 2.d0*sum0M + 2.d0*sumP0 )
               cov(3,3) = cov(3,3) + prob * ( sumPP + sumMM - 2.d0*sumPM - (sumP - sumM)**2 )
               cov(3,8) = cov(3,8) + prob * ( 1.d0/sqrt(3.d0)*( sumPP - sumMM - 2.d0*sumP0 + 2.d0*sum0M )-&
                        1.d0/sqrt(3.d0)*(sumP - sumM)*(sumP - 2.d0*sum0 + sumM) )
               cov(4,4) = cov(4,4) + prob * ( sumM + sumP + 2.d0*sumPM )
               cov(5,5) = cov(5,5) + prob * 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M - 4.d0*sumPM00 + 2.d0*sumP0 )
               cov(6,6) = cov(6,6) + prob * 0.5d0 * ( sumM + 2.d0*sum0 + sumP + 2.d0*sum0M - 4.d0*sumPM00 + 2.d0*sumP0 )
               cov(7,7) = cov(7,7) + prob * ( sumM  + sumP + 2.d0*sumPM )
               cov(8,8) = cov(8,8) + prob * ( 1.d0/3.d0*( sumPP + 4.d0*sum00 + sumMM - 4.d0*sumP0 - 4.d0*sum0M + 2.d0*sumPM )-&
                        1.d0/3.d0*(sumP - 2.d0*sum0 + sumM)**2 )

            endif

         enddo

         ! =======================================
         ! Off-diagonal part of covariance matrix
         ! =======================================
         sumJz2 = 0.0d0
         sumY2  = 0.0d0
         sumJzY = 0.0d0

         do idxi = 1, size( rhoin%eigval )

            do idxj = idxi+1, size( rhoin%eigval )

               prob  = rhoin%eigval(idxi) + rhoin%eigval(idxj)

               if( prob > 0.0_dp ) then

                  probt = rhoin%eigval(idxi) * rhoin%eigval(idxj) / prob

                  sumJz = 0.d0
                  sumY = 0.d0

                  do k=kmin, kmax

                     ! index shift
                     i = k-kmin+1

                     sumJz = sumJz + dble(M) * rhoin%eigvec(i,idxi) * rhoin%eigvec(i,idxj)
                     sumY  = sumY + 2.d0/sqrt(3.d0) * dble(3*k-N+1) * rhoin%eigvec(i,idxi) * rhoin%eigvec(i,idxj)

                  enddo

                  sumJz2 = sumJz2 + 4.d0 * probt * ( sumJz * sumJz )
                  sumY2  = sumY2 +  4.d0 * probt * ( sumY * sumY )
                  sumJzY = sumJzY + 4.d0 * probt * ( sumJz * sumY )
           
               endif

            enddo
         enddo

         cov(3,3) = cov(3,3) - sumJz2
         cov(8,8) = cov(8,8) - sumY2
         cov(3,8) = cov(3,8) - sumJzY

      endif

   end function

   ! -----------------------------------------------
   !  Optimal Fisher information for thermal state
   !     within fixed magnetization subspace
   ! -----------------------------------------------
   function Mrho_FisherQ( rhoin, N, M ) result (res)

      implicit none
      type( Mrho ), intent(in) :: rhoin 
      integer, intent(in) :: N, M
      type( fis ) :: res

      ! for lapack diagonalization
      integer, parameter :: nn = 8
      integer, parameter :: lda = 8
      integer info, lwork
      double precision a(lda, nn), w(nn)
      double precision :: work(100)

      a = dble( Mrho_Covariance( rhoin, N, M ) )

      ! ======================================
      ! Diagonalizationo of covariance matrix
      ! ======================================

      ! query the optimal workspace
      lwork = -1
      call DSYEV( 'V', 'U', nn, a, lda, w, work, lwork, info )
      lwork = int( work(1) )

      ! solve eigenvalue
      call DSYEV( 'V', 'U', nn, a, lda, w, work, lwork, info )

      res%val = w(8)/dble(N)**2
      res%dir = a(:,8)

   end function

   ! -----------------------------------------------
   !  Optimal Fisher information for thermal state
   !     with fluctuating magnetization
   ! -----------------------------------------------
   function rho_FisherQ( rhoin, N ) result (res)

      implicit none
      type( rho ), intent(in) :: rhoin 
      integer, intent(in) :: N
      type( fis ) :: res

      ! for lapack diagonalization
      integer, parameter :: nn = 8
      integer, parameter :: lda = 8
      integer info, lwork
      double precision a(lda, nn), w(nn)
      double precision :: work(100)

      a = dble( rho_Covariance( rhoin, N ) )

      ! ======================================
      ! Diagonalizationo of covariance matrix
      ! ======================================

      ! query the optimal workspace
      lwork = -1
      call DSYEV( 'V', 'U', nn, a, lda, w, work, lwork, info )
      lwork = int( work(1) )

      ! solve eigenvalue
      call DSYEV( 'V', 'U', nn, a, lda, w, work, lwork, info )

      res%val = w(8)/dble(N)**2
      res%dir = a(:,8)

   end function

   ! new covarince matrix
   function rho_Covariance( rhoin, N ) result(cov)

      implicit none
      type( rho ), intent(in) :: rhoin
      integer, intent(in) :: N
      real( kind=dp ) :: cov(8,8)

      ! local variables
      double precision sumMP, sum0P, sumM0, sumYY, sum_test, pp, qq
      integer k, kmin, kmax, i, idx, idxi, idxj, pmin, pmax, lmin, lmax, p, zmin, zmax
      integer M, Mp, Mpp
      real( kind=dp ) p1, p2, probt

      ! flush covariance matrix
      cov(:,:) = 0.0_dp

      ! loop over magnetizations
      do M = -N, N

         Mp  = M + 1  ! connecting subspace M -> Mp
         Mpp = M + 2  ! connecting subspace M -> Mpp

         if( modulo(M,2) == 0 ) then
            call EvenBounds( N, M, kmin, kmax )
            call OddBounds( N, Mp, lmin, lmax )
            call EvenBounds( N, Mpp, zmin, zmax )
         else
            call OddBounds( N, M, kmin, kmax )
            call EvenBounds( N, Mp, lmin, lmax )
            call OddBounds( N, Mpp, zmin, zmax )
         endif

         ! loop over eigenvectors
         do idxi = 1, size( rhoin%array(abs(M))%eigval )

            ! --------------------------------------------------------
            ! Loop for operator Y within the same magnetization block
            ! --------------------------------------------------------
            do idxj = idxi+1, size( rhoin%array(abs(M))%eigval )

               p1 = rhoin%array(abs(M))%eigval(idxi) * rhoin%weights(M)
               p2 = rhoin%array(abs(M))%eigval(idxj) * rhoin%weights(M)

               if( p1 + p2 > 0.d0 ) then

                  probt = ((p1 - p2)**2)/(p1 + p2)

                  sumYY = 0.d0

                  do k = kmin, kmax

                     ! index shift
                     i = k-kmin+1

                     sumYY = sumYY + dble(k) * rhoin%array(abs(M))%eigvec(i,idxi) * rhoin%array(abs(M))%eigvec(i,idxj)
   
                  enddo! k

                  sumYY = sumYY * 6.d0/sqrt(3.d0)

                  cov(8,8) = cov(8,8) + probt * sumYY**2

               endif

            enddo! idxj

            ! --------------------------------------------------------------
            ! Other operators changing magnetization
            ! --------------------------------------------------------------
            if( Mp <= N ) then

               sum_test = 0.0d0
               ! connecting M with Mp
               do idxj = 1, size( rhoin%array(abs(Mp))%eigval )

                  p1 = rhoin%array(abs(M))%eigval(idxi) * rhoin%weights(M)
                  p2 = rhoin%array(abs(Mp))%eigval(idxj) * rhoin%weights(Mp)

                  if( p1 + p2 > 0.0d0 ) then

                     probt = ((p1 - p2)**2)/(p1 + p2)
                     ! ========================
                     ! Even magnetization case
                     ! ========================
                     if( modulo(M,2) == 0 ) then

                        sum0P = 0.d0

                        pmin = max( kmin, lmin )
                        pmax = min( kmax, lmax )

                        do p = pmin, pmax
                           sum0P = sum0P + rhoin%array(abs(Mp))%eigvec(p-lmin+1,idxj) * rhoin%array(abs(M))%eigvec(p-kmin+1,idxi) *& 
                                           SQRT( dble(p + M/2 + 1) * dble(N - 2*p))
                        enddo !p

                        sumM0 = 0.d0

                        pmin = max( kmin, lmin+1 )
                        pmax = min( kmax, lmax+1 )

                        do p = pmin, pmax
                           sumM0 = sumM0 + rhoin%array(abs(Mp))%eigvec(p-lmin,idxj) * rhoin%array(abs(M))%eigvec(p-kmin+1,idxi) *&
                                           SQRT( dble(N - 2*p + 1) * dble(p - M/2) )
                        enddo !p
                  
                     ! =======================
                     ! Odd magnetization case
                     ! =======================
                     else

                        sum0P = 0.d0

                        pmin = max( kmin, lmin-1 )
                        pmax = min( kmax, lmax-1 )

                        do p = pmin, pmax
                           sum0P = sum0P + rhoin%array(abs(Mp))%eigvec(p-lmin+2,idxj) * rhoin%array(abs(M))%eigvec(p-kmin+1,idxi) *& 
                                           SQRT( dble(p + (M+1)/2 + 1) * dble(N - 2*p - 1))
                        enddo !p

                        sumM0 = 0.d0

                        pmin = max( kmin, lmin )
                        pmax = min( kmax, lmax )

                        do p = pmin, pmax
                           sumM0 = sumM0 + rhoin%array(abs(Mp))%eigvec(p-lmin+1,idxj) * rhoin%array(abs(M))%eigvec(p-kmin+1,idxi)*&
                                           SQRT( dble(N - 2*p) * dble(p - (M-1)/2) )
                        enddo !p


                     endif! if( modulo(M,2) == 0)
                  
                     sum_test = sum_test + sum0P

                     if( sum0P /= 0.0d0 ) then
                        pp = rhoin%array(abs(M))%eigval(idxi)
                        qq = rhoin%array(abs(M))%eigval(idxj)
                     endif

                     ! Jx-Jx
                     cov(1,1) = cov(1,1) + probt * abs( 1.d0/sqrt(2.d0) * (sum0P + sumM0) )**2
                     ! Jx-Qzx
                     cov(1,6) = cov(1,6) + probt * 1.d0/sqrt(2.d0) * (sum0P + sumM0) * &
                                                                  1.d0/sqrt(2.d0) * (sum0P - sumM0)
                     ! Jy-Jy
                     cov(2,2) = cov(2,2) + probt * abs( 1.d0/sqrt(2.d0) * (sum0P + sumM0) )**2
                     ! Jy-Qyz
                     cov(2,5) = cov(2,5) + probt * 1.d0/sqrt(2.d0) * (sum0P + sumM0) * &
                                                                  1.d0/sqrt(2.d0) * (sum0P - sumM0)
                     ! Qyz-Qyz
                     cov(5,5) = cov(5,5) + probt * abs( 1.d0/sqrt(2.d0) * (sum0P - sumM0) )**2
                     ! Qzx-Qzx
                     cov(6,6) = cov(6,6) + probt * abs( 1.d0/sqrt(2.d0) * (sum0P - sumM0) )**2

                  endif! (p1+p2 > 0.d0)

               enddo ! idxj

               !print *, 'M, idxi, sum, prob = ', M, idxi, sum_test, pp, qq

               ! subspace M -> Mpp
               if( Mpp <= N ) then

                  do idxj = 1, size( rhoin%array(abs(Mpp))%eigval )

                     p1 = rhoin%array(abs(M))%eigval(idxi) * rhoin%weights(M)
                     p2 = rhoin%array(abs(Mpp))%eigval(idxj) * rhoin%weights(Mpp)
                    
                     if( p1 + p2 > 0.0d0 ) then

                        probt = ((p1 - p2)**2)/(p1 + p2)
                        ! ========================
                        ! Even magnetization case
                        ! ========================
                        if( modulo(M,2) == 0 ) then

                           sumMP = 0.d0

                           pmin = max( kmin, zmin )
                           pmax = min( kmax, zmax )

                           do p = pmin, pmax
                              sumMP = sumMP + rhoin%array(abs(Mpp))%eigvec(p-zmin+1,idxj) * rhoin%array(abs(M))%eigvec(p-kmin+1,idxi)*& 
                                              SQRT( dble(p - M/2) * dble(p + M/2 + 1) )
                           enddo !p

                        ! =======================
                        ! Odd magnetization case
                        ! =======================
                        else

                           sumMP = 0.d0

                           pmin = max( kmin, zmin )
                           pmax = min( kmax, zmax )

                           do p = pmin, pmax
                              sumMP = sumMP + rhoin%array(abs(Mpp))%eigvec(p-zmin+1,idxj) * rhoin%array(abs(M))%eigvec(p-kmin+1,idxi)*& 
                                              SQRT( dble(p + (M+1)/2 + 1) * dble(p - (M-1)/2))
                           enddo !p

                        endif! (modulo(M,2) == 0)

                        ! put into covariance matrix
                        ! Qxy-Qxy
                        cov(4,4) = cov(4,4) +  probt * (sumMP)**2
                        ! Dxy-Dxy
                        cov(7,7) = cov(7,7) +  probt * (sumMP)**2

                     endif! (p1 + p2 > 0.d0)

                  enddo! idxj
               endif! (Mpp <= N)

            endif! (Mp <= N)

         enddo! idxi
      enddo! M

   end function

   ! ----------------------------------------------------
   ! fluctuations in direction 'dir' for density matrix
   ! with fixed magnetization
   ! ----------------------------------------------------
   function Mrho_Fluct( rhoin, N, M, dir ) result(fluct)

      implicit none
      type( Mrho ), intent(in) :: rhoin
      double precision, intent(in) :: dir(8)
      integer, intent(in) :: N, M
      double precision :: fluct

      integer, parameter :: nn = 8
      integer, parameter :: lda = 8
      double precision a(lda, nn), y(nn), DDOT

      a = dble( Mrho_Covariance( rhoin, N, M ) )

      call DSYMV( 'U', nn, 1.d0, a, lda, dir, 1, 0.d0, y, 1 )
      fluct = DDOT( nn, dir, 1, y, 1 )

   end function

   ! ----------------------------------------------------
   ! fluctuations in direction 'dir' for density matrix
   ! with fluctuating magnetization
   ! ----------------------------------------------------
   function rho_Fluct( rhoin, N, dir ) result(fluct)

      implicit none
      type( rho ), intent(in) :: rhoin
      double precision, intent(in) :: dir(8)
      integer, intent(in) :: N
      double precision :: fluct

      integer, parameter :: nn = 8
      integer, parameter :: lda = 8
      double precision a(lda, nn), y(nn), DDOT

      a = rho_Covariance( rhoin, N )

      call DSYMV( 'U', nn, 1.d0, a, lda, dir, 1, 0.d0, y, 1 )
      fluct = DDOT( nn, dir, 1, y, 1 )

   end function

   ! ----------------------------------------------------
   ! Fluctuations in direction 'dir' for pure state
   ! with fixed magnetization
   ! ----------------------------------------------------
   function vec_Fluct( vec, N, M, dir ) result(fluct)

      implicit none
      double precision, intent(in) :: vec(:)
      double precision, intent(in) :: dir(8)
      integer, intent(in) :: N, M
      double precision :: fluct

      integer, parameter :: nn = 8
      integer, parameter :: lda = 8
      double precision a(lda, nn), y(nn), DDOT

      a = vec_Covariance( vec, N, M )

      call DSYMV( 'U', nn, 1.d0, a, lda, dir, 1, 0.d0, y, 1 )
      Fluct = DDOT( nn, dir, 1, y, 1 )

   end function

end module
