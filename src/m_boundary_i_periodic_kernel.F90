module m_boundary_i_periodic_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_i_periodic_kernel(f,nl)
#ifdef _CUDA
      use cudafor
#endif
      use mod_dimensions, only : nx,ny,nz
      implicit none

      integer, value :: nl
      real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)

      integer :: j, k, l

#ifdef _CUDA
      j = threadIdx%y + (blockIdx%y-1)*blockDim%y-1
      k = threadIdx%z + (blockIdx%z-1)*blockDim%z-1

      if (j < 0 .or. j > ny+1) return
      if (k < 0 .or. k > nz+1) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(l,j,k) SHARED(f)
      do k = 0, nz+1
      do j = 0, ny+1
#endif

         do l= 1,nl
            f(l,0   ,j,k) = f(l,nx,j,k)
            f(l,nx+1,j,k) = f(l,1 ,j,k)
         enddo

#ifndef _CUDA
      enddo
      enddo
!$OMP END PARALLEL DO
#endif

   end subroutine
end module

