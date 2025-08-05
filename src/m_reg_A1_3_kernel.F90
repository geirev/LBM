module m_reg_A1_3_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_A1_3_kernel(A1_2, A1_3, vel, nx, ny, nz)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value :: nx, ny, nz
   real, intent(in) :: A1_2(3,3,nx,ny,nz)
   real, intent(out) :: A1_3(3,3,3,nx,ny,nz)
   real, intent(in)  :: vel(3,nx,ny,nz)
   integer :: i, j, k, p, q, r
#ifdef _CUDA
   attributes(device) :: A1_2
   attributes(device) :: A1_3
   attributes(device) :: vel
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO collapse(1) DEFAULT(NONE) PRIVATE(i, j, k, p, q, r) SHARED(vel, A1_2, A1_3, nx, ny, nz)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      do r=1,3
      do q=1,3
      do p=1,3
         A1_3(p,q,r,i,j,k)=vel(p,i,j,k)*A1_2(q,r,i,j,k) + vel(q,i,j,k)*A1_2(r,p,i,j,k) +  vel(r,i,j,k)*A1_2(p,q,i,j,k)
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
