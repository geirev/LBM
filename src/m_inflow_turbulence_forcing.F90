module m_inflow_turbulence_forcing
contains

subroutine inflow_turbulence_forcing( &
   external_forcing,rho,ampl,it,nrturb)

   use mod_dimensions, only : nx,ny,nz
   use m_inflow_turbulence_init, only : &
      iturb_pos,uu,vv,ww
#ifdef _CUDA
   use m_readinfile, only : nty,ntz
#endif
   use m_inflow_turbulence_forcing_kernel
   use m_wtime

   implicit none

   integer, intent(in) :: nrturb
   integer, intent(in) :: it
   real,    intent(in) :: ampl

   real, intent(inout) :: &
      external_forcing(3,0:nx+1,0:ny+1,0:nz+1)

   real, intent(in) :: &
      rho(0:nx+1,0:ny+1,0:nz+1)

#ifdef _CUDA
   attributes(device) :: external_forcing
   attributes(device) :: rho
#endif

   integer :: lit
   integer :: ip

#ifdef _CUDA
   integer :: tx,ty
   integer :: bx,by
#endif

   integer, parameter :: icpu=3

   call cpustart()

   ip = iturb_pos

   lit = mod(it,nrturb)
   if (lit == 0) lit=nrturb

#ifdef _CUDA
   tx=nty
   ty=ntz

   bx=(ny+tx-1)/tx
   by=(nz+ty-1)/ty
#endif

   call inflow_turbulence_forcing_kernel &
#ifdef _CUDA
      <<<dim3(bx,by,1),dim3(tx,ty,1)>>> &
#endif
      (external_forcing,rho,uu,vv,ww, &
       ampl,ip,lit)

   call cpufinish(icpu)

end subroutine inflow_turbulence_forcing

end module
