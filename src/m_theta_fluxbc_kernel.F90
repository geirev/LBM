module m_theta_fluxbc_kernel
contains

#ifdef _CUDA
attributes(global) &
#endif
subroutine theta_fluxbc_kernel(temp, heating)
   use mod_dimensions, only : nx, ny, nz
#ifdef _CUDA
   use cudafor
#endif
   implicit none

   real, intent(inout) :: temp(0:nx+1,0:ny+1,0:nz+1)
   real,         value :: heating

   integer :: i, j

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x-1)*blockDim%x
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   if (i > nx .or. j > ny) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
   do j = 1, ny
   do i = 1, nx
#endif
      ! Heat flux boundary condition at k = 0
      !
      ! dθ/dz |_{z=0} = -q0 / κ
      ! θ(0) = θ(1) + (q0/κ) * dz

      temp(i,j,0) = temp(i,j,1) + heating

#ifndef _CUDA
   end do
   end do
#endif

end subroutine theta_fluxbc_kernel

end module m_theta_fluxbc_kernel

