module m_turbines_forcing_kernel_C
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine turbines_forcing_kernel_C(turbine_df,dfeq1,dfeq2,nturbines,n)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions,  only : ny,nz
   use mod_D3Q27setup,  only : nl
   use m_turbines_init, only : ieps
   implicit none
   integer, value       :: nturbines
   integer, value       :: n
   integer, parameter   :: ii=nl*(2*ieps+1)*ny*nz
   real, intent(out)    :: turbine_df(ii,nturbines)
   real, intent(in)     :: dfeq1(ii)
   real, intent(in)     :: dfeq2(ii)
   integer i

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   if (i > ii) return
#else
!$OMP PARALLEL DO PRIVATE(i,n) SHARED(turbine_df,dfeq1,dfeq2)
      do i=1,ii
#endif
         turbine_df(i,n)=dfeq2(i)-dfeq1(i)
#ifndef _CUDA
      enddo
#endif

end subroutine
end module

