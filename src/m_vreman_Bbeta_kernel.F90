module m_vreman_Bbeta_kernel
! Eq (11) from Jacob 2018 is identical to the 33a from Feng (2021)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine vreman_Bbeta_kernel(Bbeta, beta, nx, ny, nz, nl)
   implicit none
   integer, value      :: nx, ny, nz, nl
   real, intent(out)    :: Bbeta(nx,ny,nz)
   real, intent(in)   :: beta(3,3,nx,ny,nz)
   integer :: i, j, k
#ifdef _CUDA
   attributes(device) :: Bbeta
   attributes(device) :: beta
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx) return
   if (j > ny) return
   if (k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k) SHARED(Bbeta, beta, nx, ny, nz)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
       Bbeta(i,j,k)=beta(1,1,i,j,k)*beta(2,2,i,j,k) - beta(1,2,i,j,k)**2  &
                   +beta(1,1,i,j,k)*beta(3,3,i,j,k) - beta(1,3,i,j,k)**2  &
                   +beta(2,2,i,j,k)*beta(3,3,i,j,k) - beta(2,3,i,j,k)**2
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module

