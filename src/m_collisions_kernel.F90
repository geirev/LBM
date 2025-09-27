module m_collisions_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine collisions_kernel(feq, f, tau)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions,  only : nx,ny,nz
   use mod_D3Q27setup, only : nl
   implicit none
   integer, parameter :: ntot=(nx+2)*(ny+2)*(nz+2)
   real, intent(inout) :: feq(nl*ntot)
   real, intent(in)    :: f(nl*ntot)
   real, intent(in)    :: tau(ntot)
   integer :: idx, i, j, k, it
   real :: fac

#ifdef _CUDA
   idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   if (idx > ntot*nl) return
#else
!$OMP PARALLEL DO PRIVATE(idx, i, j, k, it, fac) SHARED(f, feq, tau)
   do idx = 1, ntot*nl
#endif
     !l = mod(idx-1, nl) + 1
      i = mod( ((idx-1)/nl), nx+2 )        ! 0..nx+1
      j = mod( ((idx-1)/nl)/(nx+2), ny+2 ) ! 0..ny+1
      k = ((idx-1)/nl)/( (nx+2)*(ny+2) )   ! 0..nz+1
      it = (idx-1)/nl + 1

#ifdef _CUDA
      if (i == 0 .or. i == nx+1) return
      if (j == 0 .or. j == ny+1) return
      if (k == 0 .or. k == nz+1) return
#else
      if (i == 0 .or. i == nx+1) cycle
      if (j == 0 .or. j == ny+1) cycle
      if (k == 0 .or. k == nz+1) cycle
#endif

      fac = 1.0 - 1.0/tau(it)
      feq(idx) = feq(idx) + fac * f(idx)

#ifndef _CUDA
   end do
!$OMP END PARALLEL DO
#endif

end subroutine
end module

