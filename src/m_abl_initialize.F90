module m_abl_initialize
   real, parameter :: pottemp0    = 300.0    ! typical θ at surface [K]
contains
subroutine abl_initialize(pottemp,ir)
! Initialize potential temperature profile for ABL
! istable =  1  stable (SBL)
!            0  neutral
!           -1  unstable (convective CBL)
! pottemp returned in Kelvin
   use mod_dimensions
   use m_readinfile, only : p2l,ablheight,iablvisc,istable
   implicit none
   real, intent(out)  :: pottemp(0:nx+1,0:ny+1,0:nz+1)
   integer, intent(in) :: ir
#ifdef _CUDA
   attributes(device) :: pottemp
#endif
   real :: pottemp_h(0:nx+1,0:ny+1,0:nz+1)

! Physical parameters
   real :: dz
   real :: gamma_s      ! stable lapse rate (K/m)
   real :: gamma_n      ! neutral residual gradient
   real :: gamma_cbl    ! unstable residual gradient
   real :: gamma_top    ! gradient above BL
   real :: delta_cbl    ! inversion jump for unstable CBL
   integer :: i,j,k,is
   real :: z
   character(len=1) cstag
   if (iablvisc /=2) return

! Choose default parameters
   dz        = p2l%length
   gamma_s   = 0.015    ! 15 K/km = 0.015 K/m
   gamma_n   = 0.001    ! 1 K/km  = 0.001 K/m
   gamma_cbl = 0.001
   gamma_top = 0.005
   delta_cbl = 2.0      ! CBL inversion jump [K]


   ! Loop and fill pottemp
   do k = 0, nz+1
      z = real(k) * dz

      select case (istable)

!  STABLE ABL: θ increases strongly with height
      case (1)
         if (z <= ablheight) then
            pottemp_h(:,:,k) = pottemp0 + gamma_s * z
         else
            pottemp_h(:,:,k) = pottemp0 + gamma_s * ablheight + gamma_top * (z - ablheight)
         endif

!  NEUTRAL ABL: θ nearly constant, weak background slope
      case (0)
         pottemp_h(:,:,k) = pottemp0 + gamma_n * z

!  UNSTABLE / CONVECTIVE ABL: well-mixed layer + inversion
      case (-1)
         if (z <= ablheight) then
            pottemp_h(:,:,k) = pottemp0 + gamma_cbl*z
         else
            ! linear inversion above h
            pottemp_h(:,:,k) = pottemp0 + gamma_cbl*ablheight + delta_cbl + gamma_top*(z - ablheight)
         end if

      case default
         stop 'unvalid value for istable'
      end select
   enddo
!   pottemp_h=300.0
!   do k=0,nz+1
!   do j=0,ny+1
!   do i=0,nx+1
!      if (( (i-nx/2)**2 + (j-ny/2)**2 + (k-nz/2)**2)  < 10**2 ) then
!         pottemp_h(i,j,k)=301.0
!      endif
!   enddo
!   enddo
!   enddo

   pottemp(:,:,:)=pottemp_h(:,:,:)

   if (ir == 0) then
      is=istable
      if (istable == -1) is=2
      write(cstag,'(i1.1)')is
      open(10,file='pottemp_'//cstag//'.dat')
         do k=0,nz+1
            z = real(k) * dz
            write(10,'(i5,2f12.4)')k,z,pottemp_h(1,1,k)
         enddo
      close(10)
   endif

end subroutine
end module
