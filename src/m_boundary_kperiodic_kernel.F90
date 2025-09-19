module m_boundary_kperiodic_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_kperiodic_kernel(f,nx,ny,nz,nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value    :: nx, ny, nz, nl
   real, intent(inout) :: f(nl,nx+2,ny+2,nz+2)
   integer :: i, j, l
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x-1)*blockDim%x
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   l = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (i > nx+2) return
   if (j > ny+2) return
   if (l > nl  ) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(i,j,l) SHARED(f, nx, ny, nz, nl)
   do j=1,ny+2
   do i=1,nx+2
   do l=1,nl
#endif
      f(l,i,j,1)   =f(l,i,j,nz+1)
      f(l,i,j,nz+2)=f(l,i,j,2)
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
end subroutine
end module
