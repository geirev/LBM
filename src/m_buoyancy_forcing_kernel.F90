module m_buoyancy_forcing_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine buoyancy_forcing_kernel(external_forcing,pottemp,scaling,theta0)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions,  only : nx,ny,nz
   implicit none
   integer, parameter :: ntot=(nx+2)*(ny+2)*(nz+2)
   real, intent(inout):: external_forcing(3,ntot)
   real, intent(in)   :: pottemp(ntot)
   real, value :: scaling
   real, value :: theta0
   integer :: i,j,k,idx


#ifdef _CUDA
   idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x + 1
   if (idx > ntot) return
#else
!$OMP PARALLEL DO PRIVATE(idx, i, j, k) SHARED(external_forcing, pottemp)
   do idx = 1, ntot
#endif
      i = mod( (idx-1), nx+2 )        ! 0..nx+1
      j = mod( (idx-1)/(nx+2), ny+2 ) ! 0..ny+1
      k = (idx-1)/( (nx+2)*(ny+2) )   ! 0..nz+1

#ifdef _CUDA
      if (i == 0 .or. i == nx+1) return
      if (j == 0 .or. j == ny+1) return
      if (k == 0 .or. k == nz+1) return
#else
      if (i == 0 .or. i == nx+1) cycle
      if (j == 0 .or. j == ny+1) cycle
      if (k == 0 .or. k == nz+1) cycle
#endif

      external_forcing(3,idx) = external_forcing(3,idx) - scaling * (pottemp(idx) - theta0) / theta0

#ifndef _CUDA
   end do
!$OMP END PARALLEL DO
#endif

end subroutine
end module


