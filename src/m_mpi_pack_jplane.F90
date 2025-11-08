module m_mpi_pack_jplane
contains
#ifdef _CUDA
attributes(global) &
#endif
subroutine mpi_pack_jplane(f, j_plane, buf)
   use mod_dimensions, only : nx, ny, nz
   use mod_D3Q27setup, only : nl
   implicit none
   integer, value :: j_plane
   real :: f(nl,0:nx+1,0:ny+1,0:nz+1)   ! (l,i,j,k)
   real :: buf(nl*(nx+2)*(nz+2))       ! linearized (i,k,l); size = nl*(nx+2)*(nz+2)
#ifdef _CUDA
      attributes(device)  :: f
      attributes(device)  :: buf
#endif
   integer :: i, k, l
   integer :: idx

#ifdef _CUDA
   ! Cover i=0..nx+1, k=0..nz+1, l=1..nl  (CUDA Fortran indices are 1-based)
   i = (blockIdx%x-1)*blockDim%x + threadIdx%x - 1
   k = (blockIdx%y-1)*blockDim%y + threadIdx%y - 1
   l = (blockIdx%z-1)*blockDim%z + threadIdx%z
   if (i<0 .or. i>nx+1) return
   if (k<0 .or. k>nz+1) return
   if (l<1 .or. l>nl  ) return
#else
   do k=0,nz+1
   do i=0,nx+1
   do l=1,nl
#endif
      idx = l + nl * ( i + (nx+2) * k )
      buf(idx) = f(l,i,j_plane,k)
#ifndef _CUDA
      enddo
      enddo
      enddo
#endif
end subroutine
end module

