module m_vreman_alphamag_kernel
! Eq (11) from Jacob 2018 is identical to the 33a from Feng (2021)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine vreman_alphamag_kernel(alphamag, alpha, nx, ny, nz, nl)
   implicit none
   integer, value      :: nx, ny, nz, nl
   real, intent(in)    :: alpha(3,3,nx,ny,nz)
   real, intent(out)   :: alphamag(nx,ny,nz)
   integer :: i, j, k, l, q, p
#ifdef _CUDA
   attributes(device) :: alpha
   attributes(device) :: alphamag
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx) return
   if (j > ny) return
   if (k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, p, q) SHARED(alpha,alphamag)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      alphamag(i,j,k)=0.00001
      do q=1,3
      do p=1,3
         alphamag(i,j,k)=alphamag(i,j,k)+alpha(p,q,i,j,k)*alpha(p,q,i,j,k)
      enddo
      enddo
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module

