module m_turbine_point_forces_kernel
contains
! CUDA kernel: one thread per actuator point. Interpolates local flow, computes AoA, cl/cd, and final force vector.
#ifdef _CUDA
attributes(global) &
#endif
subroutine turbine_point_forces_kernel(points, np, rho, u, v, w, j_start, j_end, Fvec)
   use mod_dimensions, only : nx, ny, nz
   use mod_turbines,   only : point_t
   use m_nrelliftdrag
   use m_turbine_compute_bladeforce
   use m_turbine_interpolation
   implicit none

   type(point_t) :: points(:)
   integer, value :: np
   real :: rho(0:nx+1,0:ny+1,0:nz+1)
   real :: u  (0:nx+1,0:ny+1,0:nz+1)
   real :: v  (0:nx+1,0:ny+1,0:nz+1)
   real :: w  (0:nx+1,0:ny+1,0:nz+1)
   integer, value :: j_start, j_end
   real :: Fvec(3, np)

   integer :: p, j0
   real :: ux, uy, uz, dens
   real :: angattack, clift, cdrag
   real :: costheta, sintheta
   real :: utheta, phi

#ifdef _CUDA
   p = (blockIdx%x - 1) * blockDim%x + threadIdx%x
   if (p < 1 .or. p > np) return
#else
   do p=1,np
#endif

   j0 = floor(points(p)%yg)
#ifdef MPI
   if (j0 < j_start .or. j0 > j_end) then
      Fvec(1,p) = 0.0
      Fvec(2,p) = 0.0
      Fvec(3,p) = 0.0
#ifdef _CUDA
      return
#else
      cycle
#endif
   end if
#endif

   call turbine_interpolate_velocity(u, v, w, rho, &
        points(p)%xg, points(p)%yg, points(p)%zg, j_start, ux, uy, uz, dens)

   costheta = cos(points(p)%theta)
   sintheta = sin(points(p)%theta)

   utheta = points(p)%omegand * points(p)%relm - &
            uz*costheta - uy*sintheta

   phi       = atan2(ux, utheta)
   angattack = (phi*180.0/pi - points(p)%twist - points(p)%pitch)

   call nrelliftdrag(clift, cdrag, angattack, points(p)%foil)

   call turbine_compute_bladeforce(Fvec(:,p), points(p), ux, utheta, dens, clift, cdrag)
#ifndef _CUDA
   enddo
#endif
end subroutine
end module
