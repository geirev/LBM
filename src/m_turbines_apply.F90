module m_turbines_apply
contains
subroutine turbines_apply(f,F_turb,rho,u,v,w)
   use mod_dimensions,  only : nx,ny,nz
   use mod_D3Q27setup,  only : nl,cs2,cs4,cs6
   use m_readinfile   , only : ipos,nturbines,ibgk
#ifdef _CUDA
   use m_readinfile,    only : ntx,nty,ntz
#endif
   use mod_turbines, only : t_imin,t_imax, t_jmin,t_jmax, t_kmin,t_kmax
   use m_turbines_apply_kernel
   use m_wtime
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)        ! distribution
   real, intent(in)    :: F_turb(3,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    ::   u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    ::   v(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    ::   w(0:nx+1,0:ny+1,0:nz+1)

   real, allocatable :: F_turb_dev(:,:,:,:)
#ifdef _CUDA
   attributes(device) :: F_turb_dev
   attributes(device) :: f
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   
#endif
   integer n,ip

#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif
   real :: inv1cs2
   real :: inv2cs4
   real :: inv2cs6
   real :: inv6cs6
   real :: ratio
   integer, parameter :: icpu=8
   call cpustart()
   if (.not. allocated(F_turb_dev)) allocate(F_turb_dev(3,0:nx+1,0:ny+1,0:nz+1))
   F_turb_dev=F_turb

   inv1cs2 = 1.0/(cs2)
   inv2cs4 = 1.0/(2.0*cs4)
   inv2cs6 = 1.0/(2.0*cs6)
   inv6cs6 = 1.0/(6.0*cs6)
   ratio = inv6cs6 / inv2cs4
#ifdef _CUDA
   tx = ntx; bx = (t_imax - t_imin + 1 + tx - 1) / tx
   ty = nty; by = (t_jmax - t_jmin + 1 + ty - 1) / ty
   tz = ntz; bz = (t_kmax - t_kmin + 1 + tz - 1) / tz
#endif
      call turbines_apply_kernel&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          (f, rho, u, v, w, F_turb_dev,  inv1cs2, inv2cs4, ratio, ibgk,t_imin,t_imax, t_jmin,t_jmax, t_kmin,t_kmax)
   call cpufinish(icpu)
end subroutine
end module
