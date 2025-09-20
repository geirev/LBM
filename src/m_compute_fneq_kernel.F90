module m_compute_fneq_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine compute_fneq_kernel(f, feq, ntot)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: ntot
   real, intent(inout) :: f(ntot)
   real, intent(in)    :: feq(ntot)

   integer :: i

#ifdef _CUDA
   i = threadidx%x + (blockidx%x - 1) * blockdim%x
   if (i > ntot) return
#else
!$OMP PARALLEL DO  DEFAULT(none) PRIVATE(i) SHARED(f, feq, ntot)
   do i=1,ntot
#endif
      f(i) = f(i) - feq(i)
#ifndef _CUDA
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
