module m_drift_kernel
contains

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine drift_kernel(f, feq)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx, ny, nz
   use mod_D3Q27setup, only : cxs, cys, czs, nl
   implicit none

   ! f and feq flattened in (l,i), explicit dims in j,k
   real, intent(out) :: f(nl*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(in)  :: feq(nl*(nx+2), 0:ny+1, 0:nz+1)

   integer :: idx, j, k, i, l

#ifdef _CUDA
   idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x   ! flattened li index
   j   = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k   = threadIdx%z + (blockIdx%z - 1) * blockDim%z

   if (idx .le. nl .or. idx > nl*(nx+1)) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(idx,i,l,j,k) SHARED(f,feq)
   do k = 1, nz
   do j = 1, ny
   do idx = nl+1, nl*(nx+1)
#endif
      l  = mod(idx-1, nl) + 1
      i  = (idx-1)/nl
      f(idx,j,k) = feq(l + (i-cxs(l))*nl, j-cys(l), k-czs(l))
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
