module m_sphere
contains
subroutine sphere(blanking,ipos,jpos,kpos,radii)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer, intent(in)    :: ipos
   integer, intent(in)    :: jpos
   integer, intent(in)    :: kpos
   integer, intent(in)    :: radii
   integer i,j,k

   do k=1,nz
   do j=1,ny
   do i=1,nx
      if ( (i-ipos)**2 + (j-jpos)**2 + ((k-kpos)/0.5)**2 <  radii**2)  blanking(i,j,k) = .true.
   enddo
   enddo
   enddo

end subroutine
end module

