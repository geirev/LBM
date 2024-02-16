module m_channel
contains
subroutine channel(blanking,ipos,jpos,radii)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer, intent(in)    :: ipos
   integer, intent(in)    :: jpos
   integer, intent(in)    :: radii
   integer i,j

   blanking(:,1:2,:)=.true.
   blanking(:,ny-1:ny,:)=.true.

end subroutine
end module

