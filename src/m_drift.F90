module m_drift
contains
subroutine drift(f,feq)
! f enter routine in feq following collisions
! f is returned in f
   use mod_dimensions
   use mod_D3Q27setup
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime
   use m_drift_kernel
   implicit none
   real, intent(out) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)  :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer, parameter :: icpu=12
   call cpustart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   tx=ntx; bx=(nl*nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=ntz; bz=(nz+tz-1)/tz
#endif
   call drift_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f,feq)
   call cpufinish(icpu)

end subroutine
end module
