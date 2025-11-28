module m_boundary_j_periodic_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_j_periodic_kernel(f,nl)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   implicit none
   integer, value :: nl
   integer :: ntot
   real, intent(inout) :: f(ntot,ny+2,nz+2)
   integer :: i, k

   ntot=nl*(nx+2)

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x-1)*blockDim%x
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (i > ntot) return
   if (k > nz+2) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,k) SHARED(f)
   do k=1,nz+2
   do i=1,ntot
#endif
     f(i,1,k)    = f(i,ny+1,k)
     f(i,ny+2,k) = f(i,2,k)
#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
end subroutine
end module
