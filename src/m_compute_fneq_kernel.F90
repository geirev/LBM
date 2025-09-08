module m_compute_fneq_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine compute_fneq_kernel(f, feq, nx, ny, nz, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx, ny, nz, nl
   real, intent(inout) :: f(nl,nx+2,ny+2,nz+2)
   real, intent(inout) :: feq(nl,nx+2,ny+2,nz+2)



   integer :: i, j, k, l, i1, j1, k1
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
   i = threadidx%x + (blockidx%x - 1) * blockdim%x
   j = threadidx%y + (blockidx%y - 1) * blockdim%y
   k = threadidx%z + (blockidx%z - 1) * blockdim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(i, j, k, l, i1, j1, k1)&
!$OMP             & SHARED(f, feq, nx, ny, nz, nl)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      i1=i+1
      j1=j+1
      k1=k+1

      do l=1,nl
         f(l,i1,j1,k1) = f(l,i1,j1,k1) - feq(l,i1,j1,k1)
      enddo

#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
