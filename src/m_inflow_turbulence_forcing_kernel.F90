module m_inflow_turbulence_forcing_kernel
contains

#ifdef _CUDA
   attributes(global) &
#endif
subroutine inflow_turbulence_forcing_kernel( &
   external_forcing,rho,uu,vv,ww,ampl,ip,lit)

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions, only : nx,ny,nz

   implicit none

   real, intent(inout) :: &
      external_forcing(3,0:nx+1,0:ny+1,0:nz+1)

   real, intent(in) :: &
      rho(0:nx+1,0:ny+1,0:nz+1)

   real, intent(in) :: uu(ny,nz,*)
   real, intent(in) :: vv(ny,nz,*)
   real, intent(in) :: ww(ny,nz,*)

   real,    value :: ampl
   integer, value :: ip
   integer, value :: lit

   integer :: j,k

#ifdef _CUDA
   j = threadIdx%x + (blockIdx%x-1)*blockDim%x
   k = threadIdx%y + (blockIdx%y-1)*blockDim%y

   if (j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(j,k) &
!$OMP SHARED(external_forcing,rho,uu,vv,ww,ampl,ip,lit)
   do k=1,nz
   do j=1,ny
#endif

      external_forcing(1,ip,j,k) = &
         external_forcing(1,ip,j,k) &
         - rho(ip,j,k)*ampl*uu(j,k,lit)

      external_forcing(2,ip,j,k) = &
         external_forcing(2,ip,j,k) &
         - rho(ip,j,k)*ampl*vv(j,k,lit)

      external_forcing(3,ip,j,k) = &
         external_forcing(3,ip,j,k) &
         - rho(ip,j,k)*ampl*ww(j,k,lit)

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine inflow_turbulence_forcing_kernel

end module
