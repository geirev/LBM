module m_inipert_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine inipert_kernel(rho, u, v, w, rho0, stddev, udir, uvel, rho_local)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   implicit none
   real, intent(out) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(out) ::   u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(out) ::   v(0:nx+1,0:ny+1,0:nz+1)
   real, intent(out) ::   w(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) ::   rho_local(nx,ny,nz)

   real, intent(in)  :: uvel(nz)
   real, value       :: udir
   real, value       :: stddev
   real, value       :: rho0

   real, parameter   :: pi=3.1415927410125732
   integer :: i, j, k

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 1 .or. i > nx) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(u, v, w, rho, rho0, uvel, udir, stddev, pi, rho_local)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      u(i,j,k) = uvel(k) * cos(udir*pi/180.0)
      v(i,j,k) = uvel(k) * sin(udir*pi/180.0)
      w(i,j,k) = 0.0
      rho(i,j,k) = rho0 + stddev * rho_local(i,j,k)
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
end subroutine
end module

