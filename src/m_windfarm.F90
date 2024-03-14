module m_windfarm
contains
subroutine windfarm(blanking,jpos,kpos,radii)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer, intent(in)    :: jpos
   integer, intent(in)    :: kpos
   integer, intent(in)    :: radii
   integer i,j,k


   do k=1,nz
   do j=1,ny
      if ( ((j-jpos)**2 + (k-kpos)**2 ) <  radii**2) blanking(nx/4,j,k) = .true.
   enddo
   enddo
   blanking(nx/4,jpos-1:jpos+1,kpos-1:kpos+1)=.false.

   blanking(:,:,1)=.true.

end subroutine
end module

