module m_collisions
contains
subroutine collisions(f,feq,tau)
   use mod_dimensions
   implicit none
   real, intent(inout) :: f(nx,ny,nl)
   real, intent(in)    :: feq(nx,ny,nl)
   real, intent(in)    :: tau
   integer i,j,l

!$OMP PARALLEL DO PRIVATE(l,j,i) SHARED(f, feq)
      do l=1,nl
         do j=1,ny
         do i=1,nx
            f(i,j,l) = f(i,j,l) - (1.0/tau)*(f(i,j,l) - feq(i,j,l))
         enddo
         enddo
      enddo
!$OMP END PARALLEL DO
end subroutine
end module
