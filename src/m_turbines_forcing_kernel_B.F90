module m_turbines_forcing_kernel_B
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine turbines_forcing_kernel_B(u,v,w,rho,rtmp,vel,du,dv,dw,ip)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions,  only : nx,ny,nz
   use m_turbines_init, only : ieps
   implicit none
   real, intent(in)     :: rho(nx,ny,nz)
   real, intent(in)     :: u(nx,ny,nz)
   real, intent(in)     :: v(nx,ny,nz)
   real, intent(in)     :: w(nx,ny,nz)
   real, intent(out)    :: du(-ieps:ieps,ny,nz)
   real, intent(out)    :: dv(-ieps:ieps,ny,nz)
   real, intent(out)    :: dw(-ieps:ieps,ny,nz)
   real, intent(out)    :: rtmp(-ieps:ieps,ny,nz)
   real, intent(out)    :: vel(3,-ieps:ieps,ny,nz)
   integer, value       :: ip
   integer, parameter   :: ii=2*ieps+1

   integer i,j,k

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x - ieps -1
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (j > ny .or. k > nz) return
   if (i < -ieps .or. i > ieps) return
   if (ip+i < 1 .or. ip+i> nx) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(du, dv, dw, rho, ip, u, v, w, rtmp, vel)
      do k=1,nz
      do j=1,ny
      do i=-ieps,ieps
#endif
         rtmp(i,j,k) =rho(ip+i,j,k)
         vel(1,i,j,k)=u(ip+i,j,k)+du(i,j,k)
         vel(2,i,j,k)=v(ip+i,j,k)+dv(i,j,k)
         vel(3,i,j,k)=w(ip+i,j,k)+dw(i,j,k)
#ifndef _CUDA
      enddo
      enddo
      enddo
#endif

end subroutine
end module
