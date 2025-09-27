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
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use mod_D3Q27setup, only : nl
   use m_collisions_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   real, intent(in)    :: f(  nl,0:nx+1,0:ny+1,0:nz+1)  ! non-equlibrium distribution R(fneq)
   real, intent(inout) :: feq(nl,0:nx+1,0:ny+1,0:nz+1)  ! equilibrium distribution on input
   real, intent(in)    :: tau(   0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: tau
   attributes(device) :: f
   attributes(device) :: feq
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer, parameter :: ntot=(nx+2)*(ny+2)*(nz+2)
   integer, parameter :: icpu=7


   call cpustart()
#ifdef _CUDA
   tx=ntx; bx=(ntot*nl+tx-1)/tx
   ty=1; by=1
   tz=1; bz=1
#endif
   call collisions_kernel&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(feq, f, tau)
   call cpufinish(icpu)

end subroutine
end module

