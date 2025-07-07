module m_cylinder
contains
subroutine cylinder(lsolids,blanking)
   use mod_dimensions
   use m_tecfld
   logical, intent(inout) :: blanking(0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(managed) :: blanking
#endif
   logical, intent(out) :: lsolids
   integer ipos
   integer jpos
   integer radius
   integer i,j

   lsolids=.true.

   radius=5
   ipos=50
   jpos=ny/2

   do j=1,ny
   do i=1,nx
      if ( ((i-ipos)**2 + (j-jpos)**2 ) <  radius**2) blanking(i,j,0:nz) = .true.
   enddo
   enddo



end subroutine
end module
