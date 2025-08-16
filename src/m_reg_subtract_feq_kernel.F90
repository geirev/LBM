module m_reg_subtract_feq_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_subtract_feq_kernel(f, feq, nx2, ny2, nz2, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value :: nx2, ny2, nz2, nl
   real, contiguous, intent(inout) :: f(nl, nx2, ny2, nz2)
   real, contiguous, intent(in)    :: feq(nl, nx2, ny2, nz2)
   integer :: i, j, k, l
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
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(feq, f, nx2, ny2, nz2, nl)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
!! !$CUF UNROLL
!! dir$ unroll
!      do l=1,nl
!          f(l, i, j, k) = f(l, i, j, k) - feq(l, i, j, k)
!      enddo
      f( 1, i, j, k) = f( 1, i, j, k) - feq( 1, i, j, k)
      f( 2, i, j, k) = f( 2, i, j, k) - feq( 2, i, j, k)
      f( 3, i, j, k) = f( 3, i, j, k) - feq( 3, i, j, k)
      f( 4, i, j, k) = f( 4, i, j, k) - feq( 4, i, j, k)
      f( 5, i, j, k) = f( 5, i, j, k) - feq( 5, i, j, k)
      f( 6, i, j, k) = f( 6, i, j, k) - feq( 6, i, j, k)
      f( 7, i, j, k) = f( 7, i, j, k) - feq( 7, i, j, k)
      f( 8, i, j, k) = f( 8, i, j, k) - feq( 8, i, j, k)
      f( 9, i, j, k) = f( 9, i, j, k) - feq( 9, i, j, k)
      f(10, i, j, k) = f(10, i, j, k) - feq(10, i, j, k)
      f(11, i, j, k) = f(11, i, j, k) - feq(11, i, j, k)
      f(12, i, j, k) = f(12, i, j, k) - feq(12, i, j, k)
      f(13, i, j, k) = f(13, i, j, k) - feq(13, i, j, k)
      f(14, i, j, k) = f(14, i, j, k) - feq(14, i, j, k)
      f(15, i, j, k) = f(15, i, j, k) - feq(15, i, j, k)
      f(16, i, j, k) = f(16, i, j, k) - feq(16, i, j, k)
      f(17, i, j, k) = f(17, i, j, k) - feq(17, i, j, k)
      f(18, i, j, k) = f(18, i, j, k) - feq(18, i, j, k)
      f(19, i, j, k) = f(19, i, j, k) - feq(19, i, j, k)
#ifndef D3Q19
      f(20, i, j, k) = f(20, i, j, k) - feq(20, i, j, k)
      f(21, i, j, k) = f(21, i, j, k) - feq(21, i, j, k)
      f(22, i, j, k) = f(22, i, j, k) - feq(22, i, j, k)
      f(23, i, j, k) = f(23, i, j, k) - feq(23, i, j, k)
      f(24, i, j, k) = f(24, i, j, k) - feq(24, i, j, k)
      f(25, i, j, k) = f(25, i, j, k) - feq(25, i, j, k)
      f(26, i, j, k) = f(26, i, j, k) - feq(26, i, j, k)
      f(27, i, j, k) = f(27, i, j, k) - feq(27, i, j, k)
#endif
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module

