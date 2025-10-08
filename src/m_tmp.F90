module m_tmp
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine tmp(u,v,w,rho,rtmp,vel,du,dv,dw,force,iradius,ip)
#ifdef _CUDA
      use cudafor
#endif
      use mod_dimensions, only : nx,ny,nz
      use m_turbines_init, only : ieps
      implicit none

      real, intent(in)     :: rho(nx,ny,nz)
      real, intent(in)     :: u(nx,ny,nz)
      real, intent(in)     :: v(nx,ny,nz)
      real, intent(in)     :: w(nx,ny,nz)
      real, intent(in)     :: force(0:2*ieps,ny,nz,3)  ! shifted to 0:2*ieps
      real, intent(out)    :: du(-ieps:ieps,ny,nz)
      real, intent(out)    :: dv(-ieps:ieps,ny,nz)
      real, intent(out)    :: dw(-ieps:ieps,ny,nz)
      real, intent(out)    :: rtmp(-ieps:ieps,ny,nz)
      real, intent(out)    :: vel(3,-ieps:ieps,ny,nz)
      integer, value       :: iradius
      integer, value       :: ip

      integer, parameter   :: ii = 2*ieps+1
      integer :: nthreads, idx
      integer :: i,j,k,ishift

      real :: rr

#ifdef _CUDA
      nthreads = ii*ny*nz
      idx = threadIdx%x + (blockIdx%x-1)*blockDim%x
      if (idx >= nthreads) return

      i      = mod(idx, ii) - ieps
      j      = mod(idx / ii, ny) + 1
      k      = idx / (ii*ny) + 1

      if (ip+i < 1 .or. ip+i > nx) return
#else
!$OMP PARALLEL DO PRIVATE(idx,i,j,k,ishift,rr) SHARED(force,du,dv,dw,rtmp,vel,rho,u,v,w,ip)
      nthreads = ii*ny*nz
      do idx = 0, nthreads-1
         i = mod(idx, ii) - ieps
         j = mod(idx / ii, ny) + 1
         k = idx / (ii*ny) + 1
         if (ip+i < 1 .or. ip+i > nx) cycle
#endif

      ! Shift i to 0:ieps using symmetry
      ishift = abs(i)

      rr = rho(ip+i,j,k)
      du(i,j,k) = -force(ishift,j,k,1)/rr
      dv(i,j,k) = -force(ishift,j,k,2)/rr
      dw(i,j,k) = -force(ishift,j,k,3)/rr

      rtmp(i,j,k) = rr
      vel(1,i,j,k) = u(ip+i,j,k)
      vel(2,i,j,k) = v(ip+i,j,k)
      vel(3,i,j,k) = w(ip+i,j,k)

#ifndef _CUDA
      end do
#endif

   end subroutine
end module

