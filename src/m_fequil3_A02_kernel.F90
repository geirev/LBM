module m_fequil3_A02_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine fequil3_A02_kernel(rho, A0_2, vel, nx, ny, nz, inv2cs4)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value :: nx, ny, nz, nl
   real, intent(in)  :: rho(nx,ny,nz)
   real, intent(out) :: A0_2(3,3,nx,ny,nz)
   real, intent(out) :: vel(3,nx,ny,nz)
   real, value       :: inv2cs4
   integer :: i, j, k, p, q, r
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: A0_2
   attributes(device) :: vel
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k, p, q) SHARED(vel, rho, inv2cs4, A0_2, nx, ny, nz)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      do q=1,3
      do p=1,3
         A0_2(p,q,i,j,k)=rho(i,j,k)*vel(p,i,j,k)*vel(q,i,j,k)*inv2cs4
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
