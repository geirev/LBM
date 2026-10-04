module m_inflow_turbulence_forcing_kernel
contains

#ifdef _CUDA
attributes(global) &
#endif
subroutine inflow_turbulence_forcing_kernel( &
   external_forcing,rho,uu,vv,ww,iturb_i,iturb_j,iturb_face, &
   nperim,nrturb,ampl,udir,lit,ibnd,jbnd)

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions, only : nx,ny,nz

   implicit none

   real, intent(inout) :: &
      external_forcing(3,0:nx+1,0:ny+1,0:nz+1)

   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)

   integer, value :: nperim,nrturb
   integer, value :: ibnd,jbnd
   real, intent(in) :: uu(nperim,nz,0:nrturb)
   real, intent(in) :: vv(nperim,nz,0:nrturb)
   real, intent(in) :: ww(nperim,nz,0:nrturb)

   integer, intent(in) :: iturb_i(nperim)
   integer, intent(in) :: iturb_j(nperim)
   integer, intent(in) :: iturb_face(nperim)

   real,    value :: ampl,udir
   integer, value :: lit

   real, parameter :: pi=acos(-1.0)
   real, parameter :: blend_width=0.70

   integer :: p,k,i,j,iface
   real :: uxdir,uydir,x,un,activation
   real :: forcefac

#ifdef _CUDA
   p = threadIdx%x + (blockIdx%x-1)*blockDim%x
   k = threadIdx%y + (blockIdx%y-1)*blockDim%y

   if (p > nperim .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE) &
!$OMP PRIVATE(p,k,i,j,iface,uxdir,uydir,x,un,activation,forcefac) &
!$OMP SHARED(external_forcing,rho,uu,vv,ww,iturb_i,iturb_j,iturb_face, &
!$OMP        nperim,nrturb,ampl,udir,lit,ibnd,jbnd)
   do k=1,nz
   do p=1,nperim
#endif

      i=iturb_i(p)
      j=iturb_j(p)
      iface=iturb_face(p)

      uxdir=cos(udir*pi/180.0)
      uydir=sin(udir*pi/180.0)

      ! Face identifiers and inward normal wind components:
      ! 1 = y-min  : +y is inflow
      ! 2 = x-max  : -x is inflow
      ! 3 = y-max  : -y is inflow
      ! 4 = x-min  : +x is inflow
      un=0.0
      select case(iface)
      case(1)
         if (jbnd == 1) un=uydir
      case(2)
         if (ibnd == 1) un=-uxdir
      case(3)
         if (jbnd == 1) un=-uydir
      case(4)
         if (ibnd == 1) un=uxdir
      end select

      ! One-sided quintic smootherstep.  The turbulence forcing
      ! smoothly decreases to zero as the flow becomes parallel
      ! to the boundary, and is zero on an outflow boundary.
      x=min(1.0,max(0.0,un/blend_width))
      activation=x*x*x*(10.0+x*(-15.0+6.0*x))

      if (activation > 0.0) then
         ! forcings_apply uses du=-Fx/rho.  Therefore the force stored
         ! here has the opposite sign of the desired turbulent increment.
         forcefac=-rho(i,j,k)*ampl*activation

         external_forcing(1,i,j,k)=external_forcing(1,i,j,k) &
                                  +forcefac*uu(p,k,lit)
         external_forcing(2,i,j,k)=external_forcing(2,i,j,k) &
                                  +forcefac*vv(p,k,lit)
         external_forcing(3,i,j,k)=external_forcing(3,i,j,k) &
                                  +forcefac*ww(p,k,lit)
      endif

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine inflow_turbulence_forcing_kernel

end module m_inflow_turbulence_forcing_kernel
