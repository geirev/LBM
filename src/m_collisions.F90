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
   real, intent(in)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)    ! non-equlibrium distribution R(fneq)
   real, intent(inout) :: feq(nl,0:nx+1,0:ny+1,0:nz+1)  ! equilibrium distribution on input
   real, intent(in)    :: tau(nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: tau
   attributes(device) :: f
   attributes(device) :: feq
#endif
   integer i,j,k
   integer, parameter :: icpu=7
   call cpustart()


#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(f, feq, tau)
#endif
   do k=1,nz
   do j=1,ny
   do i=1,nx
      feq(:,i,j,k) =  feq(:,i,j,k) + (1.0-1.0/tau(i,j,k))*f(:,i,j,k)
   enddo
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
   call cpufinish(icpu)

end subroutine
end module

