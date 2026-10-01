module m_create_uvel_shear
contains
subroutine create_uvel_shear(filename,nz,dz,zref,z0)

   implicit none

   character(len=*), intent(in) :: filename
   integer,          intent(in) :: nz
   real,             intent(in) :: dz
   real,             intent(in) :: zref
   real,             intent(in) :: z0

   integer :: k
   real    :: z, uvel
   real    :: denom
   logical :: ex

   ! If the file already exists, leave it unchanged
   inquire(file=filename,exist=ex)
   if (ex) then
      print '(a,a)', &
         'inflow: Using existing vertical profile: ',trim(filename)
      return
   endif

   print '(a,a)', &
      'inflow: Creating vertical profile: ',trim(filename)
   print '(a,f8.3,a)', &
      '        reference height = ',zref,' m'
   print '(a,es12.4,a)', &
      '        roughness length = ',z0,' m'

   ! Neutral logarithmic ABL profile, normalized such that
   !
   !                 ln[(z+z0)/z0]
   ! uvel(z) = --------------------------
   !              ln[(zref+z0)/z0]
   !
   ! so that uvel(zref) = 1.

   denom = log((zref+z0)/z0)

   open(10,file=filename,status='new',action='write')

   do k=1,nz

      ! Fluid-node height above the halfway bounce-back wall.
      ! Wall at z=0; first fluid node at dz/2.
      z = (real(k)-0.5)*dz

      uvel = log((z+z0)/z0)/denom

      write(10,'(i6,2es20.10)') k,z,uvel

   enddo

   close(10)

end subroutine create_uvel_shear
end module
