!==============================================================
!  m_turbine_points.F90
!  Construct actuator points and compute point forces (CPU/GPU)
!==============================================================
module m_turbine_points
   use mod_turbines
   use m_turbine_extend_array
   use m_turbine_interpolation
   use m_turbine_bladeforce
#ifdef MPI
   use m_mpi_decomp_init, only : j_start, j_end, mpi_rank
#endif
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains

!--------------------------------------------------------------
!  subroutine turbine_distribute_points
!
!  PURPOSE:
!    Build the global list of actuator sample points
!    (points_global) from turbine hub locations and blade
!    geometry.
!
!  CALL:
!    call turbine_distribute_points(turbines, points_global)
!--------------------------------------------------------------
subroutine turbine_distribute_points(turbines_in, points_global)
   type(turbine_t),            intent(in)  :: turbines_in(:)
   type(point_t),   allocatable, intent(out) :: points_global(:)

   integer :: it, ib, ic
   integer :: np
   real    :: e_axis(3), e1(3), e2(3), e_rot(3)
   real    :: theta, yaw, tilt
   type(point_t) :: pt

   allocate(points_global(0))
   np = 0

   do it = 1, size(turbines_in)
      yaw  = turbines_in(it)%yaw
      tilt = turbines_in(it)%tilt

      call turbine_rotor_basis(yaw, tilt, e_axis, e1, e2)

      do ib = 1, turbines_in(it)%nblades
         theta = turbines_in(it)%theta + real(ib-1)*pi2/real(turbines_in(it)%nblades)

         do ic = 1, turbines_in(it)%nchords
            pt%iturb  = it
            pt%iblade = ib
            pt%ichord = ic

            e_rot = cos(theta)*e1 + sin(theta)*e2

            pt%xg = turbines_in(it)%xhub + turbines_in(it)%relm(ic) * e_rot(1)
            pt%yg = turbines_in(it)%yhub + turbines_in(it)%relm(ic) * e_rot(2)
            pt%zg = turbines_in(it)%zhub + turbines_in(it)%relm(ic) * e_rot(3)

            pt%yaw     = yaw
            pt%tilt    = tilt
            pt%theta   = theta
            pt%relm    = turbines_in(it)%relm(ic)
            pt%dc      = turbines_in(it)%dc(ic)
            pt%chord   = turbines_in(it)%chord(ic)
            pt%twist   = turbines_in(it)%twist(ic)
            pt%pitch   = turbines_in(it)%pitchangle
            pt%foil    = turbines_in(it)%nfoil(ic)
            pt%omegand = turbines_in(it)%omegand

            call turbine_extend_array(points_global, np+1)
            np = np + 1
            points_global(np) = pt
         end do
      end do
   end do
end subroutine turbine_distribute_points


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

      call turbine_compute_blade_force(Fvec_local(:,p), points_global(p), ux, utheta, dens, clift, cdrag)
   end do
end subroutine turbine_point_forces


!--------------------------------------------------------------
!  subroutine turbine_point_forces_gpu
!
!  PURPOSE:
!    GPU wrapper: allocate device arrays for points and Fvec,
!    launch CUDA kernel over points, then copy forces back.
!--------------------------------------------------------------
#ifdef _CUDA
subroutine turbine_point_forces_gpu(points_global, rho, u, v, w, Fvec_local, np)
   use mod_dimensions, only : nx, ny, nz
#ifdef MPI
   use m_mpi_decomp_init, only : j_start, j_end
#else
   integer, parameter :: j_start = 1, j_end = ny
#endif
   implicit none

   type(point_t), intent(in) :: points_global(:)
   real, intent(in)          :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: w  (0:nx+1,0:ny+1,0:nz+1)
   attributes(device)        :: u,v,w,rho
   integer, intent(in)       :: np
   real, intent(out)         :: Fvec_local(3, np)

   type(point_t), device, allocatable :: points_d(:)
   real,          device, allocatable :: Fvec_d(:,:)

   integer :: tpb, nblocks, istat

   if (np <= 0) then
      Fvec_local = 0.0
      return
   end if

   allocate(points_d(np))
   points_d = points_global

   allocate(Fvec_d(3,np))
   Fvec_d = 0.0

   tpb     = 128
   nblocks = (np + tpb - 1) / tpb

   call turbine_point_forces_kernel<<<nblocks, tpb>>> &
        (points_d, np, rho, u, v, w, j_start, j_end, Fvec_d)

   istat = cudaDeviceSynchronize()

   Fvec_local = Fvec_d
   !print '(a,3f13.5)','Fvec_local: ',Fvec_local(1:3,10)

   deallocate(points_d, Fvec_d)
end subroutine turbine_point_forces_gpu
#endif  /* _CUDA */


!--------------------------------------------------------------
!  subroutine turbine_point_forces_kernel
!
!  PURPOSE:
!    CUDA kernel: one thread per actuator point. Interpolates
!    local flow, computes AoA, cl/cd, and final force vector.
!--------------------------------------------------------------
#ifdef _CUDA
attributes(global) &
#endif
subroutine turbine_point_forces_kernel(points, np, rho, u, v, w, j_start, j_end, Fvec)
   use mod_dimensions, only : nx, ny, nz
   use m_nrelliftdrag
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
#else
   p = 0
#endif
   if (p < 1 .or. p > np) return

   j0 = floor(points(p)%yg)
   if (j0 < j_start .or. j0 > j_end) then
      Fvec(1,p) = 0.0
      Fvec(2,p) = 0.0
      Fvec(3,p) = 0.0
      return
   end if

   call turbine_interpolate_velocity(u, v, w, rho, &
        points(p)%xg, points(p)%yg, points(p)%zg, j_start, ux, uy, uz, dens)

   costheta = cos(points(p)%theta)
   sintheta = sin(points(p)%theta)

   utheta = points(p)%omegand * points(p)%relm - &
            uz*costheta - uy*sintheta

   phi       = atan2(ux, utheta)
   angattack = (phi*180.0/pi - points(p)%twist - points(p)%pitch)

   call nrelliftdrag(clift, cdrag, angattack, points(p)%foil)

   call turbine_compute_blade_force(Fvec(:,p), points(p), ux, utheta, dens, clift, cdrag)
end subroutine

end module m_turbine_points
