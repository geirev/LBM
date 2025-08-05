module m_vreman_beta_kernel
! Eq (11) from Jacob 2018 is identical to the 33a from Feng (2021)
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine vreman_beta_kernel(beta, alpha, nx, ny, nz, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx, ny, nz, nl
   real, intent(in)    :: alpha(3,3,nx,ny,nz)
   real, intent(out)   :: beta(3,3,nx,ny,nz)
   integer :: i, j, k, l, q, p, m
#ifdef _CUDA
   attributes(device) :: alpha
   attributes(device) :: beta
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx) return
   if (j > ny) return
   if (k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, p, q, m) SHARED(alpha, beta, nx, ny, nz)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
!! beta = del^2 * alpha' * alpha
      do q=1,3
      do p=1,3
         beta(p,q,i,j,k)=0.00001
         do m=1,3
            beta(p,q,i,j,k)=beta(p,q,i,j,k)+alpha(m,p,i,j,k)*alpha(m,q,i,j,k)
         enddo
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

