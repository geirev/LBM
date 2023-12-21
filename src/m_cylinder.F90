module m_cylinder
contains
subroutine cylinder(blanking,ipos,jpos,radii)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny)
   integer, intent(in)    :: ipos
   integer, intent(in)    :: jpos
   integer, intent(in)    :: radii
   integer i,j

   do j=1,ny
   do i=1,nx
      if ( ((i-ipos)**2 + (j-jpos)**2 ) <  radii**2) blanking(i,j) = .true.
   enddo
   enddo
!   do j=jpos-radii,jpos+radii,2
!      print '(i3,tr2,101l1)',j,blanking(ipos-2*radii:ipos+2*radii,j)
!   enddo
   open(10,file='blanking.dat')
      do j=1,ny
      do i=1,nx
         if (blanking(i,j)) write(10,'(2i4,a)')i,j,' 1.0'
      enddo
      enddo
   close(10)
end subroutine
end module 
