module m_buoyancy_forcing
contains
subroutine buoyancy_forcing(external_forcing,pottemp)
   use mod_dimensions,  only : nx,ny,nz
   use m_readinfile   , only : iablvisc
#ifdef _CUDA
   use m_readinfile,    only : ntx,nty,ntz,p2l
#endif
   use m_readinfile,    only : p2l
   use m_buoyancy_forcing_kernel
   use m_wtime
   implicit none
   integer, parameter  :: ntot=(nx+2)*(ny+2)*(nz+2)
   real, intent(inout) :: external_forcing(3,ntot)
   real, intent(in)    :: pottemp(ntot)
#ifdef _CUDA
   attributes(device) :: external_forcing
   attributes(device) :: pottemp
#endif
   real :: scaling
   real :: g
   real :: theta0
   integer n,ip

#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer, parameter :: icpu=8

   call cpustart()
   theta0=300.0
   g=9.81
   scaling=g*p2l%time**2/p2l%length

#ifdef _CUDA
   tx = ntx; bx = (ntot + tx - 1) / tx
   ty = 1; by = 1
   tz = 1; bz = 1
#endif
      call buoyancy_forcing_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(external_forcing, pottemp, scaling, theta0)
   call cpufinish(icpu)
end subroutine
end module
