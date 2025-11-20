!==============================================================
!  m_turbine_points.F90
!  Construct actuator points and compute point forces (CPU/GPU)
!==============================================================
module m_turbine_point_forces
contains


!--------------------------------------------------------------
!  subroutine turbine_point_forces
!
!  PURPOSE:
!    CPU version: for each actuator point, interpolate local
!    flow, compute angle of attack, cl/cd, and resulting
!    force vector.
!
!  CALL:
!    call turbine_point_forces(points_global, rho, u, v, w, Fvec_local, np)
!--------------------------------------------------------------
subroutine turbine_point_forces(points_global, rho, u, v, w, Fvec_local, np)
   use mod_dimensions, only : nx, ny, nz
   use m_nrelliftdrag

   use mod_turbines, only : point_t
   use m_turbine_extend_array
   use m_turbine_interpolation
   use m_turbine_compute_bladeforce
#ifdef MPI
   use m_mpi_decomp_init, only : j_start, j_end, mpi_rank
#endif
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   type(point_t), intent(in) :: points_global(:)
   real, intent(in)          :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: w  (0:nx+1,0:ny+1,0:nz+1)
   integer, intent(in)       :: np
   real,    intent(out)      :: Fvec_local(3, np)

   integer :: p, j0
   real :: ux, uy, uz, dens
   real :: angattack, clift, cdrag
   real :: costheta, sintheta
   real :: utheta, phi
   integer :: my_j_start

   Fvec_local = 0.0

#ifdef MPI
   my_j_start = j_start
#else
   my_j_start = 1
#endif

   do p = 1, np
      j0 = floor(points_global(p)%yg)
#ifdef MPI
      if (j0 < j_start .or. j0 > j_end) cycle
#endif

      call turbine_interpolate_velocity(u, v, w, rho, &
           points_global(p)%xg, points_global(p)%yg, points_global(p)%zg, &
           my_j_start, ux, uy, uz, dens)

      costheta = cos(points_global(p)%theta)
      sintheta = sin(points_global(p)%theta)

      ! Tangential relative velocity (Ω r - projection of flow)
      utheta = points_global(p)%omegand * points_global(p)%relm - &
               uz*costheta - uy*sintheta

      phi       = atan2(ux, utheta)
      angattack = (phi*180.0/pi - points_global(p)%twist - points_global(p)%pitch)

      call nrelliftdrag(clift, cdrag, angattack, points_global(p)%foil)

      call turbine_compute_bladeforce(Fvec_local(:,p), points_global(p), ux, utheta, dens, clift, cdrag)
   end do
end subroutine
end module
