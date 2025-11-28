module m_advection
contains
subroutine advection(tempout,tempin,u,v,w,tau)
! f enter routine in feq following collisions
! f is returned in f
   use mod_dimensions, only : nx, ny, nz, ntracer
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime
   use m_advection_kernel
   implicit none
   real, intent(inout) :: tempin (ntracer,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: tempout(ntracer,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: v(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: w(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: tau(0:nx+1,0:ny+1,0:nz+1)

   real            :: weights(-1:1,-1:1,-1:1)
   real, parameter :: weights_h(-1:1,-1:1,-1:1) = reshape(&
                         [ 0.125, 0.250, 0.125, 0.250, 1.000, 0.250, 0.125, 0.250, 0.125, &
                           0.250, 1.000, 0.250, 1.000,-10.000, 1.000, 0.250, 1.000, 0.250, &
                           0.125, 0.250, 0.125, 0.250, 1.000, 0.250, 0.125, 0.250, 0.125 ], shape(weights_h) )
#ifdef _CUDA
   attributes(device) :: tempin
   attributes(device) :: tempout
   attributes(device) :: u,v,w,tau
   attributes(device) :: weights
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer, parameter :: icpu=30
   integer i,j,k
   call cpustart()

   weights=weights_h

#ifdef _CUDA
   tx=ntx; bx=(ntracer*nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=ntz; bz=(nz+tz-1)/tz
#endif
   call advection_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(tempin, tempout, u, v, w, weights, tau)
   call cpufinish(icpu)

end subroutine
end module
