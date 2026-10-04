module m_wall_forcing_kernel
contains

#ifdef _CUDA
   attributes(global) &
#endif
subroutine wall_forcing_kernel(external_forcing,rho,u,v,loglaw)

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions, only : nx,ny,nz

   implicit none

   integer, parameter :: ntot=(nx+2)*(ny+2)*(nz+2)

   real, intent(inout) :: external_forcing(3,ntot)

   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v(0:nx+1,0:ny+1,0:nz+1)

   real, value :: loglaw

   real, parameter :: vonkarman=0.4

   integer :: i,j,k,idx

   real :: uh
   real :: ustar
   real :: forcefac

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x-1)*blockDim%x
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y

   if (i > nx .or. j > ny) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) &
!$OMP PRIVATE(i,j,k,idx,uh,ustar,forcefac) &
!$OMP SHARED(external_forcing,rho,u,v,loglaw)
   do j=1,ny
   do i=1,nx
#endif

      k=1

      ! Linear index corresponding to external_forcing(:,i,j,k)
      idx = 1 + i + (nx+2)*j + (nx+2)*(ny+2)*k

      ! Horizontal velocity magnitude
      uh = sqrt(u(i,j,k)**2 + v(i,j,k)**2)

      if (uh > tiny(uh)) then

         ! Friction velocity in lattice units
         ustar = vonkarman*uh/loglaw

         ! Wall shear stress distributed over the first lattice cell.
         !
         ! external_forcing is stored with the sign convention used by
         ! forcings_apply:
         !
         !        du = -Fx/rho
         !
         ! Therefore Fx,Fy below point along the velocity vector,
         ! producing a velocity increment opposite to the flow.
         !
         ! dz=1 in lattice units.

         forcefac = rho(i,j,k)*ustar*ustar/uh

         external_forcing(1,idx) = external_forcing(1,idx) &
                                  + forcefac*u(i,j,k)

         external_forcing(2,idx) = external_forcing(2,idx) &
                                  + forcefac*v(i,j,k)

      endif

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine wall_forcing_kernel

end module
