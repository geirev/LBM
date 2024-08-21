module m_channel
contains
subroutine channel(blanking,ipos,jpos,radius)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer, intent(in)    :: ipos
   integer, intent(in)    :: jpos
   integer, intent(in)    :: radius
   integer i,j

   blanking(:,1,:)=.true.
   blanking(:,ny,:)=.true.

end subroutine
end module

