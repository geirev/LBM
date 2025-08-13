module m_reg_scalef_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_scalef_kernel(f, weights, nx2, ny2, nz2, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx2, ny2, nz2, nl
   real, intent(inout) :: f(nl,1:nx2,1:ny2,1:nz2)
   real, intent(in)    :: weights(nl)
   integer :: i, j, k, l
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: weights
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 2 .or. i > nx2-1) return
   if (j < 2 .or. j > ny2-1) return
   if (k < 2 .or. k > nz2-1) return
#else
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k, l) SHARED(f, weights, nx2, ny2, nz2, nl)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
!! !$CUF UNROLL
!! dir$ unroll
!      do l=1,nl
!         f(l,i,j,k)= weights(l)*f(l,i,j,k)
!      enddo
       f(1,i,j,k)= weights(1)*f(1,i,j,k)
       f(2,i,j,k)= weights(2)*f(2,i,j,k)
       f(3,i,j,k)= weights(3)*f(3,i,j,k)
       f(4,i,j,k)= weights(4)*f(4,i,j,k)
       f(5,i,j,k)= weights(5)*f(5,i,j,k)
       f(6,i,j,k)= weights(6)*f(6,i,j,k)
       f(7,i,j,k)= weights(7)*f(7,i,j,k)
       f(8,i,j,k)= weights(8)*f(8,i,j,k)
       f(9,i,j,k)= weights(9)*f(9,i,j,k)
       f(10,i,j,k)= weights(10)*f(10,i,j,k)
       f(11,i,j,k)= weights(11)*f(11,i,j,k)
       f(12,i,j,k)= weights(12)*f(12,i,j,k)
       f(13,i,j,k)= weights(13)*f(13,i,j,k)
       f(14,i,j,k)= weights(14)*f(14,i,j,k)
       f(15,i,j,k)= weights(15)*f(15,i,j,k)
       f(16,i,j,k)= weights(16)*f(16,i,j,k)
       f(17,i,j,k)= weights(17)*f(17,i,j,k)
       f(18,i,j,k)= weights(18)*f(18,i,j,k)
       f(19,i,j,k)= weights(19)*f(19,i,j,k)
#ifndef D3Q19
       f(20,i,j,k)= weights(20)*f(20,i,j,k)
       f(21,i,j,k)= weights(21)*f(21,i,j,k)
       f(22,i,j,k)= weights(22)*f(22,i,j,k)
       f(23,i,j,k)= weights(23)*f(23,i,j,k)
       f(24,i,j,k)= weights(24)*f(24,i,j,k)
       f(25,i,j,k)= weights(25)*f(25,i,j,k)
       f(26,i,j,k)= weights(26)*f(26,i,j,k)
       f(27,i,j,k)= weights(27)*f(27,i,j,k)
#endif
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
