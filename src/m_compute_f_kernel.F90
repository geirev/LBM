module m_compute_f_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine compute_f_kernel(f, feq, ntot)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: ntot
   real, intent(inout) :: f(*)    ! f_neq or Rf_neq on input.  f on output
   real, intent(in)    :: feq(*)  ! f_eq on input
   integer :: i
#ifdef _CUDA
   i = threadidx%x + (blockidx%x - 1) * blockdim%x
   if (i > ntot) return
#else
!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(i) SHARED(f, feq, ntot)
   do i=1,ntot
#endif
      f(i) = feq(i) + f(i)
#ifndef _CUDA
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
