module m_cylinder
contains
subroutine cylinder(blanking,ipos,jpos,radius)
   use mod_dimensions
   logical, intent(inout) :: blanking(0:nx+1,0:ny+1,0:nz+1)
   integer, intent(in)    :: ipos
   integer, intent(in)    :: jpos
   integer, intent(in)    :: radius
   integer i,j

   do j=0,ny+1
   do i=0,nx+1
      if ( ((i-ipos)**2 + (j-jpos)**2 ) <  radius**2) blanking(i,j,0:30) = .true.
   enddo
   enddo

end subroutine
end module
