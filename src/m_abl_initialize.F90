module m_abl_initialize
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
   real :: pottemp0       ! surface potential temperature
   real :: gamma_s      ! stable lapse rate (K/m)
   real :: gamma_n      ! neutral residual gradient
   real :: delta_cbl    ! inversion jump for unstable CBL
   integer :: i,j,k,is
   real :: z
   character(len=1) cstag
   if (iablvisc /=2) return

! Choose default parameters
   dz        = p2l%length
   pottemp0    = 300.0    ! typical θ at surface [K]
   gamma_s   = 0.015    ! 15 K/km = 0.015 K/m
   gamma_n   = 0.001    ! 1 K/km  = 0.001 K/m
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
            pottemp_h(:,:,k) = pottemp0 + gamma_s * ablheight + 0.005 * (z - ablheight)
         endif

!  NEUTRAL ABL: θ nearly constant, weak background slope
      case (0)
         pottemp_h(:,:,k) = pottemp0 + gamma_n * z

!  UNSTABLE / CONVECTIVE ABL: well-mixed layer + inversion
      case (-1)
         if (z <= ablheight) then
            pottemp_h(:,:,k) = pottemp0                      ! well mixed
         else
            ! linear inversion above h
            pottemp_h(:,:,k) = pottemp0 + delta_cbl + 0.005*(z - ablheight)
         end if

      case default
         stop 'unvalid value for istable'
      end select
   enddo

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
