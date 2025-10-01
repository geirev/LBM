module m_boundary_iperiodic_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_iperiodic_kernel(f)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(inout) :: f(nl,nx+2,ny+2,nz+2)
   integer :: j, k, l
#ifdef _CUDA
   l = threadIdx%x + (blockIdx%x-1)*blockDim%x
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (l > nl  ) return
   if (j > ny+2) return
   if (k > nz+2) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(j,k,l) SHARED(f, nx, ny, nz, nl)
   do k=1,nz+2
   do j=1,ny+2
   do l=1,nl
#endif
      f(l,1,j,k)   =f(l,nx+1,j,k)
      f(l,nx+2,j,k)=f(l,2,j,k)
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
end subroutine
end module
