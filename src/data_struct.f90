module data_struct

public
   
   ! floating point precision
   !integer, parameter :: dp = selected_real_kind(6, 37) ! 32-bit
   !integer, parameter :: dp = selected_real_kind(15, 307) ! 64-bit
   integer, parameter :: dp = selected_real_kind(33, 4931) ! 128-bit

   ! intel 3 vector csr format
   type dcsr3
      integer :: m ! number of rows
      double precision, allocatable :: val(:)
      integer, allocatable :: col(:)
      integer, allocatable :: rowind(:)
   end type 

   type fis
      double precision :: val
      double precision :: dir(8)
   end type

   type Mrho
      double precision, allocatable :: eigval(:)
      double precision, allocatable :: eigvec(:,:)
   end type

   type rho
      type( Mrho ), allocatable :: array(:)
      real( kind = dp), allocatable :: weights(:)
   end type

contains

   subroutine EvenBounds( N, M, kmin, kmax )

      implicit none
      integer, intent(in) :: N, M
      integer, intent(out) :: kmin, kmax

      ! local variables
      real, dimension(3) :: btab

      btab = (/ 0.0, real(M)/2.0, -real(M)/2.0 /)
      kmin = ceiling( maxval( btab ) )

      btab = (/ real(N)/2.0, real(N)-real(M)/2.0, real(N)+real(M)/2.0 /)
      kmax = floor( minval( btab ) )

   end subroutine

   subroutine OddBounds( N, M, kmin, kmax )

      implicit none
      integer, intent(in) :: N, M
      integer, intent(out) :: kmin, kmax

      ! local variables
      real btab(3)

      btab = (/ 0.0, -real(M+1)/2.0, real(M-1)/2.0 /)
      kmin = ceiling( maxval( btab ) )

      btab = (/ real(N-1)/2.0, real(N)-real(M+1)/2.0, real(N)+real(M-1)/2.0 /)
      kmax = floor( minval( btab ) )

   end subroutine

end module 
