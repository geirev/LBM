module m_averaging_full_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine averaging_full_kernel(u, v, w, uave, vave, wave, uave2, vave2, wave2 )
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   implicit none
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)


   real, intent(inout) :: uave(nx,ny,nz)
   real, intent(inout) :: vave(nx,ny,nz)
   real, intent(inout) :: wave(nx,ny,nz)

   real, intent(inout) :: uave2(nx,ny,nz)
   real, intent(inout) :: vave2(nx,ny,nz)
   real, intent(inout) :: wave2(nx,ny,nz)

   integer :: i, j, k

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (k > nz) return
   if (j > ny) return
   if (i > nx) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(nx,ny,nz,uave,vave,wave,uave2,vave2,wave2,u,v,w)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      uave(i,j,k)=uave(i,j,k)+u(i,j,k)
      vave(i,j,k)=vave(i,j,k)+v(i,j,k)
      wave(i,j,k)=wave(i,j,k)+w(i,j,k)

      uave2(i,j,k)=uave2(i,j,k)+u(i,j,k)**2
      vave2(i,j,k)=vave2(i,j,k)+v(i,j,k)**2
      wave2(i,j,k)=wave2(i,j,k)+w(i,j,k)**2
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
