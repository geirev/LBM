module m_collisions_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine collisions_kernel(feq, f, tau, nx, ny, nz, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx, ny, nz, nl
   real, intent(inout) :: feq(nl,nx+2,ny+2,nz+2)
   real, intent(in)    :: f(nl,nx+2,ny+2,nz+2)
   real, intent(in)    :: tau(nx,ny,nz)
   integer :: i, j, k, l
   real fac
#ifdef _CUDA
   attributes(device) :: feq
   attributes(device) :: f
   attributes(device) :: tau
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x + 1
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y + 1
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z + 1
   if (i > nx+1) return
   if (j > ny+1) return
   if (k > nz+1) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k,fac) SHARED(f, feq, tau, nx, ny, nz, nl)
   do k=2,nz+1
   do j=2,ny+1
   do i=2,nx+1
#endif
      fac=1.0-1.0/tau(i-1,j-1,k-1)
!      do l=1,nl
!         feq(l,i,j,k) =  feq(l,i,j,k) + fac*f(l,i,j,k)
!      enddo
      feq( 1,i,j,k) =  feq( 1,i,j,k) + fac*f( 1,i,j,k)
      feq( 2,i,j,k) =  feq( 2,i,j,k) + fac*f( 2,i,j,k)
      feq( 3,i,j,k) =  feq( 3,i,j,k) + fac*f( 3,i,j,k)
      feq( 4,i,j,k) =  feq( 4,i,j,k) + fac*f( 4,i,j,k)
      feq( 5,i,j,k) =  feq( 5,i,j,k) + fac*f( 5,i,j,k)
      feq( 6,i,j,k) =  feq( 6,i,j,k) + fac*f( 6,i,j,k)
      feq( 7,i,j,k) =  feq( 7,i,j,k) + fac*f( 7,i,j,k)
      feq( 8,i,j,k) =  feq( 8,i,j,k) + fac*f( 8,i,j,k)
      feq( 9,i,j,k) =  feq( 9,i,j,k) + fac*f( 9,i,j,k)
      feq(10,i,j,k) =  feq(10,i,j,k) + fac*f(10,i,j,k)
      feq(11,i,j,k) =  feq(11,i,j,k) + fac*f(11,i,j,k)
      feq(12,i,j,k) =  feq(12,i,j,k) + fac*f(12,i,j,k)
      feq(13,i,j,k) =  feq(13,i,j,k) + fac*f(13,i,j,k)
      feq(14,i,j,k) =  feq(14,i,j,k) + fac*f(14,i,j,k)
      feq(15,i,j,k) =  feq(15,i,j,k) + fac*f(15,i,j,k)
      feq(16,i,j,k) =  feq(16,i,j,k) + fac*f(16,i,j,k)
      feq(17,i,j,k) =  feq(17,i,j,k) + fac*f(17,i,j,k)
      feq(18,i,j,k) =  feq(18,i,j,k) + fac*f(18,i,j,k)
      feq(19,i,j,k) =  feq(19,i,j,k) + fac*f(19,i,j,k)
#ifndef D3Q19
      feq(20,i,j,k) =  feq(20,i,j,k) + fac*f(20,i,j,k)
      feq(21,i,j,k) =  feq(21,i,j,k) + fac*f(21,i,j,k)
      feq(22,i,j,k) =  feq(22,i,j,k) + fac*f(22,i,j,k)
      feq(23,i,j,k) =  feq(23,i,j,k) + fac*f(23,i,j,k)
      feq(24,i,j,k) =  feq(24,i,j,k) + fac*f(24,i,j,k)
      feq(25,i,j,k) =  feq(25,i,j,k) + fac*f(25,i,j,k)
      feq(26,i,j,k) =  feq(26,i,j,k) + fac*f(26,i,j,k)
      feq(27,i,j,k) =  feq(27,i,j,k) + fac*f(27,i,j,k)
#endif
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
