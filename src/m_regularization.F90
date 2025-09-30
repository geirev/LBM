module m_regularization
contains

subroutine regularization(f, feq, u, v, w)
   use mod_dimensions
   use mod_D3Q27setup
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime

#ifdef _CUDA
    use cudafor
#endif
   use m_regularization_kernel

   implicit none
   real, intent(in)       :: u(nx,ny,nz)
   real, intent(in)       :: v(nx,ny,nz)
   real, intent(in)       :: w(nx,ny,nz)
   real, intent(inout)    :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: f
   attributes(device) :: feq
   integer :: tx, ty, tz, bx, by, bz
#endif


   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
   integer, parameter :: ntot=nl*(nx+2)*(ny+2)*(nz+2)
   integer, parameter :: icpu=5

   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#ifdef _CUDA
   tx=ntx; bx=(nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=ntz; bz=(nz+tz-1)/tz
#endif
   call regularization_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, feq, u, v, w, h2, h3, weights, inv2cs4, inv6cs6)

   call cpufinish(icpu)

end subroutine

end module
