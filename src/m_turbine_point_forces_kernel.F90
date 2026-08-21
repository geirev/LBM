module m_turbine_point_forces_kernel
contains
! CUDA kernel: one thread per actuator point. Interpolates local flow, computes AoA, cl/cd, and final force vector.
#ifdef _CUDA
attributes(global) &
#endif
subroutine turbine_point_forces_kernel(points, np, rho, u, v, w, j_start, j_end, Fvec)
   use mod_dimensions, only : nx, ny, nz
   use mod_turbines,   only : point_t,pi
   use m_liftdrag
   use m_turbine_compute_bladeforce
   use m_turbine_interpolate_velocity
   use m_turbine_rotor_basis
   implicit none

   type(point_t) :: points(:)
   integer, value :: np
   real :: rho(0:nx+1,0:ny+1,0:nz+1)
   real :: u  (0:nx+1,0:ny+1,0:nz+1)
   real :: v  (0:nx+1,0:ny+1,0:nz+1)
   real :: w  (0:nx+1,0:ny+1,0:nz+1)
   integer, value :: j_start, j_end
   real :: Fvec(3, np)

   integer :: p, ig, jg, kg
   real :: ux, uy, uz, dens
   real :: angattack, clift, cdrag
   real :: costheta, sintheta
   real :: utheta, phi
   real :: e_axis(3), e1(3), e2(3), e_tan(3)
   real :: uaxis,yaw,tilt

#ifdef _CUDA
   p = (blockIdx%x - 1) * blockDim%x + threadIdx%x
   if (p < 1 .or. p > np) return
#else
   do p=1,np
#endif

   ig = floor(points(p)%xg)
   jg = floor(points(p)%yg)
   kg = floor(points(p)%zg)

   ! Point must have a valid interpolation stencil: in i/k the upper
   ! neighbor may use the nx+1/nz+1 ghost cell, in j ownership is
   ! determined by the global lower cell against this tile's range.
   if (ig < 1 .or. ig > nx .or. &
       kg < 1 .or. kg > nz .or. &
       jg < j_start .or. jg > j_end) then
      Fvec(1,p) = 0.0
      Fvec(2,p) = 0.0
      Fvec(3,p) = 0.0
#ifdef _CUDA
      return
#else
      cycle
#endif
   end if

   call turbine_interpolate_velocity(u, v, w, rho, &
        points(p)%xg, points(p)%yg, points(p)%zg, j_start, &
        ux, uy, uz, dens)

   yaw=points(p)%yaw
   tilt=points(p)%tilt
   call turbine_rotor_basis(yaw,tilt,e_axis,e1,e2)

   costheta = cos(points(p)%theta)
   sintheta = sin(points(p)%theta)

   e_tan(:) = -sintheta*e1(:) + costheta*e2(:)

   uaxis = ux*e_axis(1) + uy*e_axis(2) + uz*e_axis(3)

   utheta = points(p)%omegand*points(p)%relm + &
            ux*e_tan(1) + uy*e_tan(2) + uz*e_tan(3)

   phi       = atan2(uaxis, utheta)
   angattack = phi*180.0/pi - points(p)%twist - points(p)%pitch

   call liftdrag(clift, cdrag, angattack, points(p)%foil)

   call turbine_compute_bladeforce(Fvec(:,p), points(p), &
                                   uaxis, utheta, dens, clift, cdrag)

   Fvec(1,p) = Fvec(1,p) * points(p)%force_scale
   Fvec(2,p) = Fvec(2,p) * points(p)%force_scale
   Fvec(3,p) = Fvec(3,p) * points(p)%force_scale

#ifndef _CUDA
   enddo
#endif
end subroutine
end module
