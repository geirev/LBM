module m_boundary_ij_corner_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_ij_corner_kernel(f)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   integer :: k, l

#ifdef _CUDA
   l = threadIdx%x + (blockIdx%x-1)*blockDim%x
   k = threadIdx%y + (blockIdx%y-1)*blockDim%y

   if (l < 1 .or. l > nl) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(l,k)
   do k = 1, nz
   do l = 1, nl
#endif

      f(l,0,0,k) = 0.5 *  (f(l,0,1,k) + f(l,1,0,k))

      f(l,0,ny+1,k) = 0.5 * (f(l,0,ny,k) + f(l,1,ny+1,k))

      f(l,nx+1,0,k) = 0.5 * (f(l,nx+1,1,k) + f(l,nx,0,k))

      f(l,nx+1,ny+1,k) = 0.5 * (f(l,nx+1,ny,k) + f(l,nx,ny+1,k))

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine boundary_ij_corner_kernel
end module m_boundary_ij_corner_kernel
