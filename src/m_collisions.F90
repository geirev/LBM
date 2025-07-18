module m_collisions
contains
subroutine collisions(f,feq,tau,it)
! returns f in feq after collisions
! NOTE:  f^coll = f - (1/tau) * (f - f^eq)
!               = f^eq + (f -f^eq) - (1/tau) * (f - f^eq)
!               = f^eq + (1-1/tau) * f^neq       # f^neq= f-f^eq
!               ~ f^eq + (1-1/tau) * R(f^neq)
   use mod_dimensions
   use m_wtime
   use m_collisions_kernel
   implicit none
   real, intent(in)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)    ! non-equlibrium distribution R(fneq)
   real, intent(inout) :: feq(nl,0:nx+1,0:ny+1,0:nz+1)  ! equilibrium distribution on input
   real, intent(in)    :: tau(nx,ny,nz)
   integer, intent(in) :: it
#ifdef _CUDA
   attributes(device) :: tau
   attributes(device) :: f
   attributes(device) :: feq
#endif
   integer i,j,k
   integer :: tx, ty, tz, bx, by, bz
   integer, parameter :: icpu=7


   call cpustart()
!!  !@cuf istat = cudaDeviceSynchronize()
!!     t0 = wallclock()
!!  #ifdef _CUDA
!!     tx=ntx; bx=(nx+2+tx-1)/tx
!!     ty=nty;  by=(ny+2+ty-1)/ty
!!     tz=ntz;  bz=(nz+2+tz-1)/tz
!!  #endif
!!     call collisions_kernel&
!!  #ifdef _CUDA
!!            &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
!!  #endif
!!            &(feq, f, tau, nx+2, ny+2, nz+2, nl)
!!  
!!  
!!  

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

!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(31)=walltimelocal(31)+t1-t0

   call cpufinish(icpu)

   if (it==999) then
      do j=31,31
         print '(a24,i3,g13.5)','collisions:',j,walltimelocal(j)
      enddo
   endif

end subroutine
end module

