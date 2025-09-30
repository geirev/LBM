module m_macrovars
! Computes density as a sum over particles with different velocities
contains
subroutine macrovars(rho,u,v,w,f)
   use mod_dimensions
   use mod_D3Q27setup
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime
   use m_macrovars_kernel
   implicit none
   real,    intent(in)  :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(out) :: rho(nx,ny,nz)
   real,    intent(out) :: u(nx,ny,nz)
   real,    intent(out) :: v(nx,ny,nz)
   real,    intent(out) :: w(nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer, parameter :: icpu=13

   call cpustart()
#ifdef _CUDA
   tx=ntx; bx=(nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=ntz; bz=(nz+tz-1)/tz
#endif
   call macrovars_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, rho, u, v, w)

   call cpufinish(icpu)

end subroutine
end module
