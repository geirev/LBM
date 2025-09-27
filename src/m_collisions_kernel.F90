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
   real, intent(inout) :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    ::   f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    ::    tau(0:nx+1,0:ny+1,0:nz+1)
   integer :: i, j, k , l
   real fac
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx) return
   if (j > ny) return
   if (k > nz) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k,fac) SHARED(f, feq, tau, nx, ny, nz, nl)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      fac=1.0-1.0/tau(i,j,k)
      do l=1,nl
         feq(l,i,j,k) =  feq(l,i,j,k) + fac*f(l,i,j,k)
      enddo

#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
