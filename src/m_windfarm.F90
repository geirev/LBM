module m_windfarm
contains
subroutine windfarm(blanking,ipos,jpos,kpos,radii)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer, intent(in)    :: ipos    ! i index of turbine
   integer, intent(in)    :: jpos    ! j index of turbine
   integer, intent(in)    :: kpos    ! k index of turbine
   integer, intent(in)    :: radii   ! turbine radius
   integer i,j,k


! Adding a trubine
   do k=1,nz
   do j=1,ny
      if ( ((j-jpos)**2 + (k-kpos)**2 ) <  radii**2) blanking(ipos,j,k) = .true.
   enddo
   enddo
   blanking(ipos,jpos-2:jpos+2,kpos-2:kpos+2)=.false. ! adding a hole at turbine center

! No-slip at ground
   blanking(:,:,1)=.true.

end subroutine
end module

