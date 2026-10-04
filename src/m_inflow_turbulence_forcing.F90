module m_inflow_turbulence_forcing
contains

subroutine inflow_turbulence_forcing( &
   external_forcing,rho,ampl,udir,it,nrturb)

   use mod_dimensions, only : nx,ny,nz
   use m_inflow_turbulence_init, only : &
      nturb_perim,iturb_i,iturb_j,iturb_face,uu,vv,ww
#ifdef _CUDA
   use m_readinfile, only : nty,ntz
#endif
   use m_readinfile, only : ibnd,jbnd
   use m_inflow_turbulence_forcing_kernel
   use m_wtime

   implicit none

   integer, intent(in) :: nrturb,it
   real,    intent(in) :: ampl,udir

   real, intent(inout) :: &
      external_forcing(3,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)

#ifdef _CUDA
   attributes(device) :: external_forcing,rho
#endif

   integer :: lit

#ifdef _CUDA
   integer :: tx,ty,bx,by
#endif

   integer, parameter :: icpu=3

   call cpustart()
   if (ibnd /= 1 .and. jbnd /= 1) return

   lit=mod(it,nrturb)
   if (lit == 0) lit=nrturb

#ifdef _CUDA
   tx=nty
   ty=ntz
   bx=(nturb_perim+tx-1)/tx
   by=(nz+ty-1)/ty
#endif

   call inflow_turbulence_forcing_kernel &
#ifdef _CUDA
      <<<dim3(bx,by,1),dim3(tx,ty,1)>>> &
#endif
      (external_forcing,rho,uu,vv,ww, &
       iturb_i,iturb_j,iturb_face,nturb_perim,nrturb, &
       ampl,udir,lit,ibnd,jbnd)

   call cpufinish(icpu)

end subroutine inflow_turbulence_forcing

end module m_inflow_turbulence_forcing
