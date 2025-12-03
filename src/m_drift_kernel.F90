module m_drift_kernel
contains

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine drift_kernel(f, feq)
#ifdef _CUDA
      use cudafor
#endif
      use mod_dimensions,  only : nx, ny, nz
      use mod_D3Q27setup, only : cxs, cys, czs, nl
      implicit none

      real, intent(inout) :: f  (nl,0:nx+1,0:ny+1,0:nz+1)
      real, intent(inout) :: feq(nl,0:nx+1,0:ny+1,0:nz+1)

      integer :: idx, i, l, j, k

#ifdef _CUDA
      ! Flattened (l,i) index: idx = 1..nl*nx
      idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x
      j   = threadIdx%y + (blockIdx%y - 1) * blockDim%y
      k   = threadIdx%z + (blockIdx%z - 1) * blockDim%z

      if (idx < 1 .or. idx > nl*nx) return
      if (j   < 1 .or. j   > ny   ) return
      if (k   < 1 .or. k   > nz   ) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(idx,i,l,j,k) SHARED(f,feq)
      do k = 1, nz
      do j = 1, ny
      do idx = 1, nl*nx
#endif
         ! Map flat idx back to (l,i):
         ! idx = (i-1)*nl + l  ->  l = mod(idx-1,nl)+1, i = (idx-1)/nl + 1
         l = mod(idx-1, nl) + 1
         i = (idx-1) / nl + 1     ! i = 1..nx

         ! Standard streaming:
         ! f(l,i,j,k) comes from feq at (i - cx, j - cy, k - cz)
         f(l,i,j,k) = feq(l, i - cxs(l), j - cys(l), k - czs(l))

#ifndef _CUDA
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
#endif

   end subroutine
end module

