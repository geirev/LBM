module m_collisions
contains
subroutine collisions(f,feq,tau)
   use mod_dimensions
   implicit none
   real, intent(inout) :: f(nx,ny,nz,nl)
   real, intent(in)    :: feq(nx,ny,nz,nl)
   real, intent(in)    :: tau
   integer l

!$OMP PARALLEL DO PRIVATE(l) SHARED(f, feq, tau)
      do l=1,nl
         f(:,:,:,l) = f(:,:,:,l) - (1.0/tau)*(f(:,:,:,l) - feq(:,:,:,l))
      enddo
!$OMP END PARALLEL DO
end subroutine
end module
