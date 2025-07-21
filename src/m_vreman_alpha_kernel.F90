module m_vreman_alpha_kernel
! Eq (11) from Jacob 2018 is identical to the 33a from Feng (2021)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine vreman_alpha_kernel(f, H2, alpha, nx2, ny2, nz2, nl)
   implicit none
   integer, value      :: nx2, ny2, nz2, nl
   real, intent(in)    :: f(nl,nx2,ny2,nz2)
   real, intent(in)    :: H2(3,3,nl)
   real, intent(out)   :: alpha(3,3,nx2-2,ny2-2,nz2-2)
   integer :: i, j, k, l, q, p
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: H2
   attributes(device) :: alpha
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 2 .or. i > nx2-1) return
   if (j < 2 .or. j > ny2-1) return
   if (k < 2 .or. k > nz2-1) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q) SHARED(f, H2, alpha, nl, nx2, ny2, nz2)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
      do q=1,3
      do p=1,3
         alpha(p,q,i-1,j-1,k-1) =  H2(p,q,1)*f(1,i,j,k)
      enddo
      enddo

      do l=2,nl
         do q=1,3
         do p=1,3
            alpha(p,q,i-1,j-1,k-1) = alpha(p,q,i-1,j-1,k-1) + H2(p,q,l)*f(l,i,j,k)
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
