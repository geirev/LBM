module m_collisions_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine collisions_kernel(feq, f, tau, ntot, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: ntot, nl
   real, intent(inout) :: feq(nl*ntot)
   real, intent(in)    :: f(nl*ntot)
   real, intent(in)    :: tau(ntot)
   integer :: idx, i, l
   real :: fac

#ifdef _CUDA
   idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   if (idx > ntot*nl) return
#else
!$OMP PARALLEL DO PRIVATE(idx, i, l, fac) SHARED(f, feq, tau, ntot, nl)
   do idx = 1, ntot*nl
#endif
      ! recover (l,i) from flattened index
      l = mod(idx-1, nl) + 1
      i = (idx-1)/nl + 1

      fac = 1.0 - 1.0 / tau(i)
      feq(idx) = feq(idx) + fac * f(idx)

#ifndef _CUDA
   end do
!$OMP END PARALLEL DO
#endif

end subroutine
end module

