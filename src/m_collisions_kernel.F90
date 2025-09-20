module m_collisions_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine collisions_kernel(feq, f, tau, ntot, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: ntot, nl
   real, intent(inout) :: feq(nl,ntot)
   real, intent(in)    :: f(nl,ntot)
   real, intent(in)    :: tau(ntot)
   integer :: i,l
   real fac
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   if (i > ntot) return
#else
!$OMP PARALLEL DO PRIVATE(i,j,k,fac) SHARED(f, feq, tau, ntot, nl)
   do i=1,ntot
#endif
      fac=1.0-1.0/tau(i)
      do l=1,nl
         feq(l,i) =  feq(l,i) + fac*f(l,i)
      enddo

#ifndef _CUDA
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
