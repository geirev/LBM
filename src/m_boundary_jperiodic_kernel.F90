module m_boundary_jperiodic_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_jperiodic_kernel(f)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(inout) :: f(nl,nx+2,ny+2,nz+2)
   integer :: i, k, l
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x-1)*blockDim%x
   l = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (i > nx+2) return
   if (l > nl  ) return
   if (k > nz+2) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,k,l) SHARED(f, nx, ny, nz, nl)
   do k=1,nz+2
   do i=1,nx+2
   do l=1,nl
#endif
     f(l,i,1,k)    = f(l,i,ny+1,k)
     f(l,i,ny+2,k) = f(l,i,2,k)
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
end subroutine
end module
