module m_atmospheric_grid

   implicit none
   private

   public :: atmospheric_grid_t
   public :: atmospheric_grid_allocate
   public :: atmospheric_interp
   public :: atmospheric_lbm_boundaries

   !---------------------------------------------------------------
   ! Atmospheric model grid.
   !
   ! x(ia), y(ja), z(ka) are physical coordinates.
   !
   ! u(ia,ja,ka), v(ia,ja,ka), w(ia,ja,ka) are the atmospheric
   ! velocity components in physical units.
   !---------------------------------------------------------------
   type :: atmospheric_grid_t

      integer :: nxa = 0
      integer :: nya = 0
      integer :: nza = 0

      real, allocatable :: x(:)
      real, allocatable :: y(:)
      real, allocatable :: z(:)

      real, allocatable :: u(:,:,:)
      real, allocatable :: v(:,:,:)
      real, allocatable :: w(:,:,:)

   end type atmospheric_grid_t

contains


!=======================================================================
! Allocate atmospheric grid
!=======================================================================

subroutine atmospheric_grid_allocate(atm,nxa,nya,nza)

   implicit none

   type(atmospheric_grid_t), intent(inout) :: atm
   integer, intent(in) :: nxa,nya,nza

   atm%nxa = nxa
   atm%nya = nya
   atm%nza = nza

   allocate(atm%x(nxa))
   allocate(atm%y(nya))
   allocate(atm%z(nza))

   allocate(atm%u(nxa,nya,nza))
   allocate(atm%v(nxa,nya,nza))
   allocate(atm%w(nxa,nya,nza))

   atm%x = 0.0
   atm%y = 0.0
   atm%z = 0.0

   atm%u = 0.0
   atm%v = 0.0
   atm%w = 0.0

end subroutine atmospheric_grid_allocate


!=======================================================================
! Interpolate atmospheric velocity at physical coordinate (xp,yp,zp).
!
! Trilinear interpolation is used.
!
! Coordinates outside the atmospheric domain are clamped to the
! closest atmospheric grid interval. Normally the atmospheric grid
! should completely enclose the LBM domain.
!=======================================================================

subroutine atmospheric_interp(atm,xp,yp,zp,up,vp,wp)

   implicit none

   type(atmospheric_grid_t), intent(in) :: atm

   real, intent(in)  :: xp,yp,zp
   real, intent(out) :: up,vp,wp

   integer :: i0,j0,k0

   real :: ax,ay,az

   real :: w000,w100,w010,w110
   real :: w001,w101,w011,w111


   ! Find surrounding atmospheric cell and local coordinates

   call find_interval(atm%x,atm%nxa,xp,i0,ax)
   call find_interval(atm%y,atm%nya,yp,j0,ay)

   if (atm%nza > 1) then

      call find_interval(atm%z,atm%nza,zp,k0,az)

   else

      k0 = 1
      az = 0.0

   endif


   ! Trilinear interpolation weights

   w000 = (1.0-ax)*(1.0-ay)*(1.0-az)
   w100 =       ax *(1.0-ay)*(1.0-az)
   w010 = (1.0-ax)*      ay *(1.0-az)
   w110 =       ax *      ay *(1.0-az)

   if (atm%nza > 1) then

      w001 = (1.0-ax)*(1.0-ay)*az
      w101 =       ax *(1.0-ay)*az
      w011 = (1.0-ax)*      ay *az
      w111 =       ax *      ay *az


      !---------------------------------------------------------------
      ! u component
      !---------------------------------------------------------------

      up = w000*atm%u(i0  ,j0  ,k0  ) + &
           w100*atm%u(i0+1,j0  ,k0  ) + &
           w010*atm%u(i0  ,j0+1,k0  ) + &
           w110*atm%u(i0+1,j0+1,k0  ) + &
           w001*atm%u(i0  ,j0  ,k0+1) + &
           w101*atm%u(i0+1,j0  ,k0+1) + &
           w011*atm%u(i0  ,j0+1,k0+1) + &
           w111*atm%u(i0+1,j0+1,k0+1)


      !---------------------------------------------------------------
      ! v component
      !---------------------------------------------------------------

      vp = w000*atm%v(i0  ,j0  ,k0  ) + &
           w100*atm%v(i0+1,j0  ,k0  ) + &
           w010*atm%v(i0  ,j0+1,k0  ) + &
           w110*atm%v(i0+1,j0+1,k0  ) + &
           w001*atm%v(i0  ,j0  ,k0+1) + &
           w101*atm%v(i0+1,j0  ,k0+1) + &
           w011*atm%v(i0  ,j0+1,k0+1) + &
           w111*atm%v(i0+1,j0+1,k0+1)


      !---------------------------------------------------------------
      ! w component
      !---------------------------------------------------------------

      wp = w000*atm%w(i0  ,j0  ,k0  ) + &
           w100*atm%w(i0+1,j0  ,k0  ) + &
           w010*atm%w(i0  ,j0+1,k0  ) + &
           w110*atm%w(i0+1,j0+1,k0  ) + &
           w001*atm%w(i0  ,j0  ,k0+1) + &
           w101*atm%w(i0+1,j0  ,k0+1) + &
           w011*atm%w(i0  ,j0+1,k0+1) + &
           w111*atm%w(i0+1,j0+1,k0+1)

   else

      !---------------------------------------------------------------
      ! Bilinear interpolation for a single atmospheric z-level
      !---------------------------------------------------------------

      up = w000*atm%u(i0  ,j0  ,1) + &
           w100*atm%u(i0+1,j0  ,1) + &
           w010*atm%u(i0  ,j0+1,1) + &
           w110*atm%u(i0+1,j0+1,1)

      vp = w000*atm%v(i0  ,j0  ,1) + &
           w100*atm%v(i0+1,j0  ,1) + &
           w010*atm%v(i0  ,j0+1,1) + &
           w110*atm%v(i0+1,j0+1,1)

      wp = w000*atm%w(i0  ,j0  ,1) + &
           w100*atm%w(i0+1,j0  ,1) + &
           w010*atm%w(i0  ,j0+1,1) + &
           w110*atm%w(i0+1,j0+1,1)

   endif

end subroutine atmospheric_interp


!=======================================================================
! Locate q in monotonically increasing coordinate array x.
!
! On return:
!
!       x(i0) <= q <= x(i0+1)
!
! and
!
!       q = (1-a)*x(i0) + a*x(i0+1)
!
! with 0 <= a <= 1.
!
! Values outside the coordinate range are clamped.
!=======================================================================

subroutine find_interval(x,n,q,i0,a)

   implicit none

   integer, intent(in) :: n

   real, intent(in) :: x(n)
   real, intent(in) :: q

   integer, intent(out) :: i0
   real,    intent(out) :: a

   integer :: i


   if (n < 2) then

      i0 = 1
      a  = 0.0

      return

   endif


   ! Lower edge

   if (q <= x(1)) then

      i0 = 1
      a  = 0.0

      return

   endif


   ! Upper edge

   if (q >= x(n)) then

      i0 = n-1
      a  = 1.0

      return

   endif


   ! Interior

   do i=1,n-1

      if (q >= x(i) .and. q <= x(i+1)) then

         i0 = i
         a  = (q-x(i))/(x(i+1)-x(i))

         return

      endif

   enddo


   ! Should never reach here

   i0 = n-1
   a  = 1.0

end subroutine find_interval


!=======================================================================
! Interpolate atmospheric velocity onto all four horizontal LBM
! boundaries.
!
! LBM physical coordinates:
!
!    x(i) = (i-1)*dx
!    y(j) = (j-1)*dx
!    z(k) = (k-1)*dx
!
! Boundary arrays:
!
!    xmin: x = 0
!    xmax: x = (nx-1)*dx
!
!       dimensions (ny,nz)
!
!    ymin: y = 0
!    ymax: y = (ny-1)*dx
!
!       dimensions (nx,nz)
!
! The coordinates are those of the PHYSICAL LBM boundary nodes,
! not the ghost nodes.
!=======================================================================

subroutine atmospheric_lbm_boundaries(atm,nx,ny,nz,dx, &
                                      u_xmin,v_xmin,w_xmin, &
                                      u_xmax,v_xmax,w_xmax, &
                                      u_ymin,v_ymin,w_ymin, &
                                      u_ymax,v_ymax,w_ymax)

   implicit none

   type(atmospheric_grid_t), intent(in) :: atm

   integer, intent(in) :: nx,ny,nz
   real,    intent(in) :: dx

   real, intent(out) :: u_xmin(ny,nz)
   real, intent(out) :: v_xmin(ny,nz)
   real, intent(out) :: w_xmin(ny,nz)

   real, intent(out) :: u_xmax(ny,nz)
   real, intent(out) :: v_xmax(ny,nz)
   real, intent(out) :: w_xmax(ny,nz)

   real, intent(out) :: u_ymin(nx,nz)
   real, intent(out) :: v_ymin(nx,nz)
   real, intent(out) :: w_ymin(nx,nz)

   real, intent(out) :: u_ymax(nx,nz)
   real, intent(out) :: v_ymax(nx,nz)
   real, intent(out) :: w_ymax(nx,nz)

   integer :: i,j,k

   real :: x,y,z

   real :: xmin,xmax
   real :: ymin,ymax


   xmin = 0.0
   xmax = real(nx-1)*dx

   ymin = 0.0
   ymax = real(ny-1)*dx


   !---------------------------------------------------------------
   ! x boundaries
   !---------------------------------------------------------------

   do k=1,nz

      z = real(k-1)*dx

      do j=1,ny

         y = real(j-1)*dx

         call atmospheric_interp(atm,xmin,y,z, &
                                 u_xmin(j,k),   &
                                 v_xmin(j,k),   &
                                 w_xmin(j,k))

         call atmospheric_interp(atm,xmax,y,z, &
                                 u_xmax(j,k),   &
                                 v_xmax(j,k),   &
                                 w_xmax(j,k))

      enddo

   enddo


   !---------------------------------------------------------------
   ! y boundaries
   !---------------------------------------------------------------

   do k=1,nz

      z = real(k-1)*dx

      do i=1,nx

         x = real(i-1)*dx

         call atmospheric_interp(atm,x,ymin,z, &
                                 u_ymin(i,k),   &
                                 v_ymin(i,k),   &
                                 w_ymin(i,k))

         call atmospheric_interp(atm,x,ymax,z, &
                                 u_ymax(i,k),   &
                                 v_ymax(i,k),   &
                                 w_ymax(i,k))

      enddo

   enddo

end subroutine atmospheric_lbm_boundaries


end module m_atmospheric_grid
