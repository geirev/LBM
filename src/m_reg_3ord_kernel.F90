module m_reg_3ord_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_3ord_kernel(f, A1_3, H3, nx, ny, nz, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value    :: nx, ny, nz, nl
   real, intent(out) :: f(nl,nx+2,ny+2,nz+2)
   real, intent(in)  :: A1_3(3,3,3,nx,ny,nz)
   real, intent(in)  :: H3(3,3,3,nl)
   integer :: i, j, k, l, p, q, r
#ifdef _CUDA
   attributes(device) :: A1_3
   attributes(device) :: H3
   attributes(device) :: f
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, r) SHARED(f, H3, A1_3, nx, ny, nz, nl)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      do l=2,nl
         do r=1,3
         do q=1,3
         do p=1,3
            f(l,i+1,j+1,k+1)=f(l,i+1,j+1,k+1) + H3(p,q,r,l)*A1_3(p,q,r,i,j,k)
         enddo
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
