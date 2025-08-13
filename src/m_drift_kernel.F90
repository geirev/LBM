module m_drift_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine drift_kernel(f,feq, nx2, ny2, nz2, nl, cxs, cys, czs)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value    :: nx2, ny2, nz2, nl
   real, intent(out) :: f(nl,nx2,ny2,nz2)
   real, intent(in)  :: feq(nl,nx2,ny2,nz2)
   integer, intent(in) :: cxs(nl)
   integer, intent(in) :: cys(nl)
   integer, intent(in) :: czs(nl)
   integer :: i, j, k, l
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
   attributes(device) :: cxs,cys,czs
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 2 .or. i > nx2-1) return
   if (j < 2 .or. j > ny2-1) return
   if (k < 2 .or. k > nz2-1) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k,l) SHARED(f, feq, nx2, ny2, nz2, nl)
   do k=2,nz2-1
   do j=2,ny2-1
   do i=2,nx2-1
#endif
!! !$CUF UNROLL
!! dir$ unroll
!      do l=1,nl
!         f(l,i,j,k) = feq(l,i-cxs(l),j-cys(l),k-czs(l))
!      enddo
       f( 1,i,j,k) = feq( 1,i  ,j  ,k  )
       f( 2,i,j,k) = feq( 2,i-1,j  ,k  )
       f( 3,i,j,k) = feq( 3,i+1,j  ,k  )
       f( 4,i,j,k) = feq( 4,i  ,j-1,k  )
       f( 5,i,j,k) = feq( 5,i  ,j+1,k  )
       f( 6,i,j,k) = feq( 6,i  ,j  ,k+1)
       f( 7,i,j,k) = feq( 7,i  ,j  ,k-1)
       f( 8,i,j,k) = feq( 8,i-1,j-1,k  )
       f( 9,i,j,k) = feq( 9,i+1,j+1,k  )
       f(10,i,j,k) = feq(10,i-1,j+1,k  )
       f(11,i,j,k) = feq(11,i+1,j-1,k  )
       f(12,i,j,k) = feq(12,i+1,j  ,k+1)
       f(13,i,j,k) = feq(13,i-1,j  ,k-1)
       f(14,i,j,k) = feq(14,i  ,j-1,k-1)
       f(15,i,j,k) = feq(15,i  ,j+1,k+1)
       f(16,i,j,k) = feq(16,i+1,j  ,k-1)
       f(17,i,j,k) = feq(17,i-1,j  ,k+1)
       f(18,i,j,k) = feq(18,i  ,j+1,k-1)
       f(19,i,j,k) = feq(19,i  ,j-1,k+1)
       f(20,i,j,k) = feq(20,i+1,j-1,k-1)
       f(21,i,j,k) = feq(21,i-1,j+1,k+1)
       f(22,i,j,k) = feq(22,i+1,j+1,k+1)
       f(23,i,j,k) = feq(23,i-1,j-1,k-1)
       f(24,i,j,k) = feq(24,i-1,j-1,k+1)
       f(25,i,j,k) = feq(25,i+1,j+1,k-1)
       f(26,i,j,k) = feq(26,i+1,j-1,k+1)
       f(27,i,j,k) = feq(27,i-1,j+1,k-1)
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
