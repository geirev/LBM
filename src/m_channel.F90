module m_channel
contains
subroutine channel(blanking,ipos,jpos,radii)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer, intent(in)    :: ipos
   integer, intent(in)    :: jpos
   integer, intent(in)    :: radii
   integer i,j

   do j=1,ny
   do i=1,nx
      if ( ((i-ipos)**2 + (j-jpos)**2 ) <  radii**2) blanking(i,j,:) = .true.
   enddo
   enddo

   blanking(:,1,:)=.true.
   blanking(:,ny,:)=.true.

end subroutine
end module

