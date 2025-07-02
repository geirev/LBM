module m_sphere
contains
subroutine sphere(lsolids,blanking,ipos,jpos,kpos,radius)
   use mod_dimensions
   logical, intent(inout) :: blanking(0:nx+1,0:ny+1,0:nz+1)
   integer, intent(in)    :: ipos
   integer, intent(in)    :: jpos
   integer, intent(in)    :: kpos
   integer, intent(in)    :: radius
   integer i,j,k
   logical, intent(inout) :: lsolids
   lsolids=.true.

   do k=0,nz+1
   do j=0,ny+1
   do i=0,nx+1
      if ( (i-ipos)**2 + (j-jpos)**2 + (k-kpos)**2 <  radius**2)  blanking(i,j,k) = .true.
   enddo
   enddo
   enddo

end subroutine
end module

