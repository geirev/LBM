module m_compute_f
contains

subroutine compute_f(f, feq)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime

#ifdef _CUDA
    use cudafor
#endif
   use m_compute_f_kernel

   implicit none
   real, intent(inout)    :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
   integer :: tx, ty, tz, bx, by, bz
#endif

   integer, parameter :: icpu=5

   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing non-equilibrium distribution defined in \citet{fen21a} between Eqs (32) and (33)
#ifdef _CUDA
   tx=ntx; bx=(nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=ntz; bz=(nz+tz-1)/tz
#endif
   call compute_f_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, feq, nx, ny, nz, nl)


   call cpufinish(icpu)

end subroutine

end module
