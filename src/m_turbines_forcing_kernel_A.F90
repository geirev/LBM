module m_turbines_forcing_kernel_A
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine turbines_forcing_kernel_A(u,v,w,rho,rtmp,vel,du,dv,dw,force,iradius,ip)
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
   real, intent(in)     :: force(0:ieps,ny,nz,3)
   real, intent(out)    :: du(-ieps:ieps,ny,nz)
   real, intent(out)    :: dv(-ieps:ieps,ny,nz)
   real, intent(out)    :: dw(-ieps:ieps,ny,nz)
   real, intent(out)    :: rtmp(-ieps:ieps,ny,nz)
   real, intent(out)    :: vel(3,-ieps:ieps,ny,nz)
   integer, value       :: iradius
   integer, value       :: ip
   integer, parameter   :: ii=2*ieps+1

   real rr
   integer i,j,k

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x - ieps -1
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (j > ny .or. k > nz) return
   if (i < -ieps .or. i > ieps) return
   if (ip+i < 1 .or. ip+i> nx) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(force, du, dv, dw, rho, iradius, u, v, w, ip, rtmp, vel)
      do k=1,nz
      do j=1,ny
      do i=-ieps,ieps
#endif
         rr=rho(ip+i,j,k)
         du(i,j,k)=-force(abs(i),j,k,1)/rr
         dv(i,j,k)=-force(abs(i),j,k,2)/rr
         dw(i,j,k)=-force(abs(i),j,k,3)/rr

         rtmp(i,j,k) =rho(ip+i,j,k)
         vel(1,i,j,k)=u(ip+i,j,k)
         vel(2,i,j,k)=v(ip+i,j,k)
         vel(3,i,j,k)=w(ip+i,j,k)
#ifndef _CUDA
      enddo
      enddo
      enddo
#endif

end subroutine
end module
