module m_collisions_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine collisions_kernel(feq, f, tau, nx2, ny2, nz2, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx2, ny2, nz2, nl
   real, intent(inout) :: feq(nl,nx2,ny2,nz2)
   real, intent(in)    :: f(nl,nx2,ny2,nz2)
   real, intent(in)    :: tau(nx2-2,ny2-2,nz2-2)
   integer :: i, j, k, l
#ifdef _CUDA
   attributes(device) :: feq
   attributes(device) :: f
   attributes(device) :: tau
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 2 .or. i > nx2-1) return
   if (j < 2 .or. j > ny2-1) return
   if (k < 2 .or. k > nz2-1) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(f, feq, tau, nx2, ny2, nz2, nl)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
         feq(:,i,j,k) =  feq(:,i,j,k) + (1.0-1.0/tau(i-1,j-1,k-1))*f(:,i,j,k)
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
