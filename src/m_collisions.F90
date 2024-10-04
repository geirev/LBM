module m_collisions
contains
subroutine collisions(f,feq,tau)
! returns f in feq after collisions
! NOTE:  f^coll = f - (1/tau) * (f - f^eq)
!               = f^eq + (f -f^eq) - (1/tau) * (f - f^eq)
!               = f^eq + (1-1/tau) * f^neq       # f^neq= f-f^eq
!               ~ f^eq + (1-1/tau) * R(f^neq)
   use mod_dimensions
   use m_wtime
   implicit none
   real, intent(in)    :: f(0:nx+1,0:ny+1,0:nz+1,nl)    ! non-equlibrium distribution R(fneq)
   real, intent(inout) :: feq(0:nx+1,0:ny+1,0:nz+1,nl)  ! equilibrium distribution on input
   real, intent(in)    :: tau(nx,ny,nz)
   integer l,i,j,k
   integer, parameter :: icpu=6
   call cpustart()


!$OMP PARALLEL DO PRIVATE(l,i,j,k) SHARED(f, feq, tau)
      do l=1,nl
         do k=1,nz
         do j=1,ny
         do i=1,nx
            feq(i,j,k,l) =  feq(i,j,k,l) + (1.0-1.0/tau(i,j,k))*f(i,j,k,l)
         enddo
         enddo
         enddo
      enddo
!$OMP END PARALLEL DO
   call cpufinish(icpu)

end subroutine
end module

