module m_compute_fneq
contains

subroutine compute_fneq(f, feq)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime

#ifdef _CUDA
    use cudafor
#endif
   use m_compute_fneq_kernel

   implicit none
   real, intent(in)    :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
   integer :: tx, ty, tz, bx, by, bz
#endif

   integer, parameter :: ntot=nl*(nx+2)*(ny+2)*(nz+2)
   integer, parameter :: icpu=5

   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing non-equilibrium distribution defined in \citet{fen21a} between Eqs (32) and (33)
#ifdef _CUDA
   tx=ntx; bx=(ntot+tx-1)/tx
   ty=1; by=1
   tz=1; bz=1
#endif
   call compute_fneq_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, feq, ntot)


   call cpufinish(icpu)

end subroutine

end module
