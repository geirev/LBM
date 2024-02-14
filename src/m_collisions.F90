module m_collisions
contains
subroutine collisions(f,feq,tau)
! returns f in feq after collisions
   use mod_dimensions
   implicit none
   real, intent(in)    :: f(0:nx+1,ny,nz,nl)
   real, intent(inout) :: feq(0:nx+1,ny,nz,nl)
   real, intent(in)    :: tau
   integer l

!$OMP PARALLEL DO PRIVATE(l) SHARED(f, feq, tau)
      do l=1,nl
         feq(:,:,:,l) =  f(:,:,:,l) - (1.0/tau)*(f(:,:,:,l) - feq(:,:,:,l))
      enddo
!$OMP END PARALLEL DO

end subroutine
end module
