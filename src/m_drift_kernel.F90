module m_drift_kernel
contains

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine drift_kernel(f, feq)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx, ny, nz
   use mod_D3Q27setup, only : cxs, cys, czs, nl
   implicit none

   ! f and feq flattened in (l,i), explicit dims in j,k
   real, intent(out) :: f(nl*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(in)  :: feq(nl*(nx+2), 0:ny+1, 0:nz+1)

   integer :: idx, j, k, i, l

#ifdef _CUDA
   idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x   ! flattened li index
   j   = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k   = threadIdx%z + (blockIdx%z - 1) * blockDim%z

   if (idx .le. nl .or. idx > nl*(nx+1)) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(idx,i,l,j,k) SHARED(f,feq)
   do k = 1, nz
   do j = 1, ny
   do idx = nl+1, nl*(nx+1)
#endif
      l  = mod(idx-1, nl) + 1
      i  = (idx-1)/nl
      f(idx,j,k) = feq(l + (i-cxs(l))*nl, j-cys(l), k-czs(l))
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module






!!  module m_drift_kernel
!!  contains
!!  #ifdef _CUDA
!!     attributes(global)&
!!  #endif
!!     subroutine drift_kernel(f,feq)
!!  #ifdef _CUDA
!!     use cudafor
!!  #endif
!!     use mod_dimensions, only : nx, ny, nz
!!     use mod_D3Q27setup, only : cxs, cys, czs, nl
!!     implicit none
!!     real, intent(out) ::   f(nl,0:nx+1,0:ny+1,0:nz+1)
!!     real, intent(in)  :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
!!     integer :: i, j, k
!!  #ifdef _CUDA
!!     i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
!!     j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
!!     k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
!!     if (i > nx) return
!!     if (j > ny) return
!!     if (k > nz) return
!!  #else
!!  !$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,k) SHARED(f, feq)
!!     do k=1,nz
!!     do j=1,ny
!!     do i=1,nx
!!  #endif
!!  ! !$CUF UNROLL
!!  !      do l=1,nl
!!  !         f(l,i,j,k) = feq(l,i-cxs(l),j-cys(l),k-czs(l))
!!  !      enddo
!!         f( 1,i,j,k) = feq( 1,i  ,j  ,k  )
!!         f( 2,i,j,k) = feq( 2,i-1,j  ,k  )
!!         f( 3,i,j,k) = feq( 3,i+1,j  ,k  )
!!         f( 4,i,j,k) = feq( 4,i  ,j-1,k  )
!!         f( 5,i,j,k) = feq( 5,i  ,j+1,k  )
!!         f( 6,i,j,k) = feq( 6,i  ,j  ,k+1)
!!         f( 7,i,j,k) = feq( 7,i  ,j  ,k-1)
!!         f( 8,i,j,k) = feq( 8,i-1,j-1,k  )
!!         f( 9,i,j,k) = feq( 9,i+1,j+1,k  )
!!         f(10,i,j,k) = feq(10,i-1,j+1,k  )
!!         f(11,i,j,k) = feq(11,i+1,j-1,k  )
!!         f(12,i,j,k) = feq(12,i+1,j  ,k+1)
!!         f(13,i,j,k) = feq(13,i-1,j  ,k-1)
!!         f(14,i,j,k) = feq(14,i  ,j-1,k-1)
!!         f(15,i,j,k) = feq(15,i  ,j+1,k+1)
!!         f(16,i,j,k) = feq(16,i+1,j  ,k-1)
!!         f(17,i,j,k) = feq(17,i-1,j  ,k+1)
!!         f(18,i,j,k) = feq(18,i  ,j+1,k-1)
!!         f(19,i,j,k) = feq(19,i  ,j-1,k+1)
!!  #ifndef D3Q19
!!         f(20,i,j,k) = feq(20,i+1,j-1,k-1)
!!         f(21,i,j,k) = feq(21,i-1,j+1,k+1)
!!         f(22,i,j,k) = feq(22,i+1,j+1,k+1)
!!         f(23,i,j,k) = feq(23,i-1,j-1,k-1)
!!         f(24,i,j,k) = feq(24,i-1,j-1,k+1)
!!         f(25,i,j,k) = feq(25,i+1,j+1,k-1)
!!         f(26,i,j,k) = feq(26,i+1,j-1,k+1)
!!         f(27,i,j,k) = feq(27,i-1,j+1,k-1)
!!  #endif
!!  #ifndef _CUDA
!!     enddo
!!     enddo
!!     enddo
!!  !$OMP END PARALLEL DO
!!  #endif
!!  
!!  end subroutine
!!  end module
