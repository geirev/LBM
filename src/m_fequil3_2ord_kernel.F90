! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
module m_fequil3_2ord_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine fequil3_2ord_kernel(feq, H2, A0_2, nx2, ny2, nz2, nl)
   implicit none
   integer, value      :: nx2, ny2, nz2, nl
   real, intent(inout) :: feq(nl,nx2,ny2,nz2)
   real, intent(in)    :: H2(3,3,nl)
   real, intent(in)    :: A0_2(3,3,nx2-2,ny2-2,nz2-2)
   integer :: i, j, k, l, q, p
#ifdef _CUDA
   attributes(device) :: feq
   attributes(constant) :: H2
   attributes(device) :: A0_2
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 2 .or. i > nx2-1) return
   if (j < 2 .or. j > ny2-1) return
   if (k < 2 .or. k > nz2-1) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q) SHARED(feq, H2, A0_2, nx2, ny2, nz2, nl)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
      do l=1,nl
         do q=1,3
         do p=1,3
            feq(l,i,j,k)=feq(l,i,j,k) + H2(p,q,l)*A0_2(p,q,i-1,j-1,k-1)
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
