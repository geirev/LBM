module m_cylinder
contains
subroutine cylinder(blanking,ipos,jpos,radius)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer, intent(in)    :: ipos
   integer, intent(in)    :: jpos
   integer, intent(in)    :: radius
   integer i,j

   do j=1,ny
   do i=1,nx
      if ( ((i-ipos)**2 + (j-jpos)**2 ) <  radius**2) blanking(i,j,1:30) = .true.
   enddo
   enddo

end subroutine
end module
