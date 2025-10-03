module m_boundary_k_periodic_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_k_periodic_kernel(f)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl
   implicit none
   integer, parameter     :: ntot=nl*(nx+2)
   real, intent(inout) :: f(ntot,ny+2,nz+2)
   integer :: i, j
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x-1)*blockDim%x
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   if (i > ntot) return
   if (j > ny+2) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j) SHARED(f)
   do j=1,ny+2
   do i=1,ntot
#endif
      f(i,j,1)   =f(i,j,nz+1)
      f(i,j,nz+2)=f(i,j,2)
#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
end subroutine
end module
