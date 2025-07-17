module m_reg_cp_vel_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_cp_vel_kernel(vel, u, v, w, nx, ny, nz)
   implicit none
   integer, value :: nx, ny, nz
   real, intent(in) :: u(nx, ny, nz)
   real, intent(in) :: v(nx, ny, nz)
   real, intent(in) :: w(nx, ny, nz)
   real, intent(out):: vel(3,nx,ny,nz)
   integer :: i, j, k
#ifdef _CUDA
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: vel
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(vel, u, v, w, nx, ny, nz)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      vel(1,i,j,k)=u(i,j,k)
      vel(2,i,j,k)=v(i,j,k)
      vel(3,i,j,k)=w(i,j,k)
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module

