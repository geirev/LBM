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
!$OMP PARALLEL DO PRIVATE(i,fac) SHARED(f, feq, tau, ntot, nl)
   do i=1,ntot
#endif
      fac=1.0-1.0/tau(i)
!      do l=1,nl
!         feq(l,i) =  feq(l,i) + fac*f(l,i)
!      enddo
      feq(1,i)  =  feq(1,i)  + fac*f(1,i)
      feq(2,i)  =  feq(2,i)  + fac*f(2,i)
      feq(3,i)  =  feq(3,i)  + fac*f(3,i)
      feq(4,i)  =  feq(4,i)  + fac*f(4,i)
      feq(5,i)  =  feq(5,i)  + fac*f(5,i)
      feq(6,i)  =  feq(6,i)  + fac*f(6,i)
      feq(7,i)  =  feq(7,i)  + fac*f(7,i)
      feq(8,i)  =  feq(8,i)  + fac*f(8,i)
      feq(9,i)  =  feq(9,i)  + fac*f(9,i)
      feq(10,i) =  feq(10,i) + fac*f(10,i)
      feq(11,i) =  feq(11,i) + fac*f(11,i)
      feq(12,i) =  feq(12,i) + fac*f(12,i)
      feq(13,i) =  feq(13,i) + fac*f(13,i)
      feq(14,i) =  feq(14,i) + fac*f(14,i)
      feq(15,i) =  feq(15,i) + fac*f(15,i)
      feq(16,i) =  feq(16,i) + fac*f(16,i)
      feq(17,i) =  feq(17,i) + fac*f(17,i)
      feq(18,i) =  feq(18,i) + fac*f(18,i)
      feq(19,i) =  feq(19,i) + fac*f(19,i)
#ifndef D3Q19
      feq(20,i) =  feq(20,i) + fac*f(20,i)
      feq(21,i) =  feq(21,i) + fac*f(21,i)
      feq(22,i) =  feq(22,i) + fac*f(22,i)
      feq(23,i) =  feq(23,i) + fac*f(23,i)
      feq(24,i) =  feq(24,i) + fac*f(24,i)
      feq(25,i) =  feq(25,i) + fac*f(25,i)
      feq(26,i) =  feq(26,i) + fac*f(26,i)
      feq(27,i) =  feq(27,i) + fac*f(27,i)
#endif

#ifndef _CUDA
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
