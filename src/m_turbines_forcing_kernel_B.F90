module m_turbines_forcing_kernel_B
! Assumes vel was set in kernel_A and not changed after
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine turbines_forcing_kernel_B(vel,du,dv,dw)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions,  only : nx,ny,nz
   use m_turbines_init, only : ieps
   implicit none
   integer, parameter :: ntot=(2*ieps+1)*ny*nz
   real, intent(out)    :: du(ntot)
   real, intent(out)    :: dv(ntot)
   real, intent(out)    :: dw(ntot)
   real, intent(inout)    :: vel(3,ntot)

   integer i

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   if (i > ntot) return
#else
!$OMP PARALLEL DO PRIVATE(i) SHARED(du, dv, dw, vel)
      do i=1,ntot
#endif
         vel(1,i)=vel(1,i)+du(i)
         vel(2,i)=vel(2,i)+dv(i)
         vel(3,i)=vel(3,i)+dw(i)
#ifndef _CUDA
      enddo
#endif

end subroutine
end module
