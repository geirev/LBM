module m_reg_scalef_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_scalef_kernel(f, weights, nx2, ny2, nz2, nl)
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
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k, l) SHARED(f, weights, nx, ny, nz, nl)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
      do l=1,nl
         f(l,i,j,k)= weights(l)*f(l,i,j,k)
         !f(l,i+1,j+1,k+1)= weights(l)*f(l,i+1,j+1,k+1)
      enddo
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
