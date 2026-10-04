module m_create_uvel_shear
contains
subroutine create_uvel_shear(filename,nz,dz,zref,z0,uvel_shear)

   implicit none

   character(len=*), intent(in)  :: filename
   integer,          intent(in)  :: nz
   real,             intent(in)  :: dz
   real,             intent(in)  :: zref
   real,             intent(in)  :: z0
   real,             intent(out) :: uvel_shear(nz)

   integer :: k
   integer :: iu
   real    :: z
   real    :: denom

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

   open(newunit=iu,file=filename,status='replace',action='write')

   ! Tecplot header
   write(iu,'(a)') 'TITLE = "Vertical inflow wind profile"'
   write(iu,'(a)') 'VARIABLES = "K", "Z", "UVEL"'
   write(iu,'(a,i0,a)') &
      'ZONE T="Inflow profile", I=',nz,', F=POINT'

   do k=1,nz

      ! Fluid-node height above the halfway bounce-back wall.
      ! Wall at z=0; first fluid node at dz/2.
      z = (real(k)-0.5)*dz

      uvel_shear(k) = log((z+z0)/z0)/denom

      write(iu,'(i6,2es20.10)') k,z,uvel_shear(k)

   enddo

   close(iu)
end subroutine create_uvel_shear
end module
