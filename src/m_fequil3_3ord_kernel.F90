! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
module m_fequil3_3ord_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine fequil3_3ord_kernel(feq, H3, A0_3, nx, ny, nz, nl)
   implicit none
   integer, value      :: nx, ny, nz, nl
   real, intent(inout) :: feq(nl,nx+2,ny+2,nz+2)
   real, intent(in)    :: H3(3,3,3,nl)
   real, intent(in)    :: A0_3(3,3,3,nx,ny,nz)
   integer :: i, j, k, l, q, p, r
#ifdef _CUDA
   attributes(device) :: feq
   attributes(device) :: H3
   attributes(device) :: A0_3
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q) SHARED(feq, H2, A0_2)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      do l=2,nl
         do r=1,3
         do q=1,3
         do p=1,3
            feq(l,i+1,j+1,k+1)=feq(l,i+1,j+1,k+1) + H3(p,q,r,l)*A0_3(p,q,r,i,j,k)
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
