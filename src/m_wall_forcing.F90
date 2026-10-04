module m_wall_forcing
contains

subroutine wall_forcing(external_forcing,rho,u,v)

   use mod_dimensions, only : nx,ny,nz
   use m_readinfile,   only : p2l,z0
#ifdef _CUDA
   use m_readinfile,   only : ntx,nty
#endif
   use m_wall_forcing_kernel
   use m_wtime

   implicit none

   integer, parameter :: ntot=(nx+2)*(ny+2)*(nz+2)

   real, intent(inout) :: external_forcing(3,ntot)
   real, intent(in)    :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: v(0:nx+1,0:ny+1,0:nz+1)

#ifdef _CUDA
   attributes(device) :: external_forcing
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
#endif

   real :: z0_lattice
   real :: loglaw

#ifdef _CUDA
   integer :: tx,ty,bx,by
#endif

   integer, parameter :: icpu=9

   call cpustart()

   ! Roughness length in lattice units
   z0_lattice = z0/p2l%length

   ! First fluid node is located at z=0.5 in lattice units
   loglaw = log((0.5 + z0_lattice)/z0_lattice)

#ifdef _CUDA
   tx = ntx
   ty = nty

   bx = (nx + tx - 1)/tx
   by = (ny + ty - 1)/ty
#endif

   call wall_forcing_kernel &
#ifdef _CUDA
      <<<dim3(bx,by,1),dim3(tx,ty,1)>>> &
#endif
      (external_forcing,rho,u,v,loglaw)

   call cpufinish(icpu)

end subroutine wall_forcing

end module
