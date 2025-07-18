module m_reg_subtract_feq_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_subtract_feq_kernel(f, feq, nx2, ny2, nz2, nl)
   implicit none
   integer, value :: nx2, ny2, nz2, nl
   real, intent(inout) :: f(nl, nx2, ny2, nz2)
   real, intent(in)    :: feq(nl, nx2, ny2, nz2)
   integer :: i, j, k
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 2 .or. i > nx2-1) return
   if (j < 2 .or. j > ny2-1) return
   if (k < 2 .or. k > nz2-1) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(feq, f, nx, ny, nz, nl)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
        f(:, i, j, k) = f(:, i, j, k) - feq(:, i, j, k)
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module

