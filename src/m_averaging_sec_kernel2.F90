module m_averaging_sec_kernel2
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine averaging_sec_kernel2(u, v, w, uxave, vxave, wxave, uxave2, vxave2, wxave2, kpos )
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   implicit none
   integer, value      :: kpos
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)


   real, intent(inout) :: uxave(nx,ny)
   real, intent(inout) :: vxave(nx,ny)
   real, intent(inout) :: wxave(nx,ny)

   real, intent(inout) :: uxave2(nx,ny)
   real, intent(inout) :: vxave2(nx,ny)
   real, intent(inout) :: wxave2(nx,ny)

   integer ::  i,j
#ifdef _CUDA
   integer ::  k
#endif

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k=1
   if (j > ny) return
   if (i > nx) return
#else
!$OMP PARALLEL DO PRIVATE(i,j) SHARED(nx, ny, u, v, w, uxave, vxave, wxave, uxave2, vxave2, wxave2, kpos )
   do j=1,ny
   do i=1,nx
#endif
      uxave(i,j)=uxave(i,j) + u(i,j,kpos)
      vxave(i,j)=vxave(i,j) + v(i,j,kpos)
      wxave(i,j)=wxave(i,j) + w(i,j,kpos)

      uxave2(i,j)=uxave2(i,j) + u(i,j,kpos)**2
      vxave2(i,j)=vxave2(i,j) + v(i,j,kpos)**2
      wxave2(i,j)=wxave2(i,j) + w(i,j,kpos)**2
#ifndef _CUDA
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
