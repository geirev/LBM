module m_reg_subtract_feq_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_subtract_feq_kernel(f, feq, nx, ny, nz, nl)
   implicit none
   integer, value :: nx, ny, nz, nl ! nx,ny,nz has values+2
   real, intent(inout) :: f(nl, nx+2, ny+2, nz+2)
   real, intent(in)    :: feq(nl, nx+2, ny+2, nz+2)
   integer :: i, j, k, l
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(feq, f, nx, ny, nz, nl)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      do l = 1, nl
        f(l, i+1, j+1, k+1) = f(l, i+1, j+1, k+1) - feq(l, i+1, j+1, k+1)
      enddo
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif


!!   #ifdef _CUDA
!!   !$cuf kernel do(2) <<<*,*>>>
!!   #else
!!   !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(feq, f)
!!   #endif
!!         do k=1,nz
!!            do j=1,ny
!!               do i=1,nx
!!                  do l=1,nl
!!                     f(l,i,j,k)=f(l,i,j,k)-feq(l,i,j,k)
!!                  enddo
!!               enddo
!!            enddo
!!         enddo

end subroutine
end module

