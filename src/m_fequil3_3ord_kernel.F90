! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
module m_fequil3_3ord_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine fequil3_3ord_kernel(feq, H3, A0_3, nx2, ny2, nz2, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx2, ny2, nz2, nl
   real, intent(inout) :: feq(nl,nx2,ny2,nz2)
   real, intent(in)    :: H3(3,3,3,nl)
   real, intent(in)    :: A0_3(3,3,3,nx2-2,ny2-2,nz2-2)
   integer :: i, j, k, l, q, p, r
#ifdef _CUDA
   attributes(device) :: feq
   attributes(device) :: H3
   attributes(device) :: A0_3
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 2 .or. i > nx2-1) return
   if (j < 2 .or. j > ny2-1) return
   if (k < 2 .or. k > nz2-1) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q) SHARED(feq, nx2, ny2, nz2, nl, H3, A0_3)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
      do l=2,nl
         do r=1,3
         do q=1,3
         do p=1,3
            feq(l,i,j,k)=feq(l,i,j,k) + H3(p,q,r,l)*A0_3(p,q,r,i-1,j-1,k-1)
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
