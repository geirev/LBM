module m_reg_scaleA1_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_scaleA1_kernel(A1_2, A1_3, nx, ny, nz, inv2cs4, inv6cs6)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx, ny, nz
   real, intent(inout) :: A1_2(3,3,nx,ny,nz)
   real, intent(inout) :: A1_3(3,3,3,nx,ny,nz)
   real, value         :: inv2cs4,inv6cs6
   integer :: i, j, k, p, q, r
#ifdef _CUDA
   attributes(device) :: A1_2
   attributes(device) :: A1_3
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO collapse(1) DEFAULT(NONE) PRIVATE(i, j, k, p, q, r) SHARED(A1_2, A1_3, inv2cs4, inv6cs6, nx, ny, nz)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif

      do q=1,3
      do p=1,3
         A1_2(p,q,i,j,k) = A1_2(p,q,i,j,k)*inv2cs4
      enddo
      enddo

      do r=1,3
      do q=1,3
      do p=1,3
         A1_3(p,q,r,i,j,k) = A1_3(p,q,r,i,j,k)*inv6cs6
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
