module m_buoyancy_forcing_apply
contains
subroutine buoyancy_forcing_apply(f,buoyancy_force,rho,u,v,w)
   use mod_dimensions,  only : nx,ny,nz
   use mod_D3Q27setup,  only : nl,cs2,cs4,cs6
   use m_readinfile   , only : iablvisc,ibgk
   use m_turbines_apply_kernel
#ifdef _CUDA
   use m_readinfile,    only : ntx,nty,ntz
#endif
   use m_wtime
   implicit none
   real, intent(inout)  :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)     :: buoyancy_force(3,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)     :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)     :: u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)     :: v(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)     :: w(0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device)   :: f,buoyancy_force,rho,u,v,w
#endif
   real :: inv1cs2
   real :: inv2cs4
   real :: inv2cs6
   real :: inv6cs6
   real :: ratio
   integer, parameter :: icpu=29
   integer t_imin, t_imax, t_jmin,t_jmax, t_kmin,t_kmax
#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif

   call cpustart()
   inv1cs2 = 1.0/(cs2)
   inv2cs4 = 1.0/(2.0*cs4)
   inv2cs6 = 1.0/(2.0*cs6)
   inv6cs6 = 1.0/(6.0*cs6)
   ratio = inv6cs6 / inv2cs4

#ifdef _CUDA
   tx = ntx; bx = (nx + tx - 1) / tx
   ty = nty; by = (ny + ty - 1) / ty
   tz = ntz; bz = (nz + tz - 1) / tz
#endif
      call turbines_apply_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(f, rho, u, v, w, buoyancy_force, inv1cs2, inv2cs4, ratio, ibgk, t_imin, t_imax, t_jmin,t_jmax, t_kmin,t_kmax)
   call cpufinish(icpu)
   end subroutine
end module


