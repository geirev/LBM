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
      real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)

      integer :: i, k, l

#ifdef _CUDA
      i = threadIdx%x + (blockIdx%x-1)*blockDim%x-1
      k = threadIdx%z + (blockIdx%z-1)*blockDim%z-1

      if (i < 0 .or. i > nx+1) return
      if (k < 0 .or. k > nz+1) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,k,l) SHARED(f)
      do k = 0, nz+1
      do i = 0, nx+1
#endif

         do l=1,nl
            f(l,i,0   ,k) = f(l,i,ny,k)
            f(l,i,ny+1,k) = f(l,i,1 ,k)
         enddo

#ifndef _CUDA
      enddo
      enddo
!$OMP END PARALLEL DO
#endif

   end subroutine
end module

