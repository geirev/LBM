module m_cylinder
contains
subroutine cylinder(blanking)
   use mod_dimensions, only : nx, nyg, nz
   implicit none
   logical, intent(inout) :: blanking(0:nx+1,0:nyg+1,0:nz+1)
   integer :: i,j
   integer, parameter :: radius = 7, ipos = 50, jpos = nyg/2

   blanking = .false.

   do j=1,nyg
      do i=1,nx
         if ((i-ipos)**2 + (j-jpos)**2 < radius**2) blanking(i,j,:) = .true.
      end do
   end do
end subroutine
end module
