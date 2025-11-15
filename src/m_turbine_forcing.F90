!==============================================================
!  m_turbine_forcing.F90
!  High-level driver to compute turbine forcing on each tile
!==============================================================
module m_turbine_forcing
contains
!--------------------------------------------------------------
!  subroutine turbine_forcing
!
!  PURPOSE:
!    High-level driver called from the main LBM time loop.
!    For each step:
!      1) Update rotor azimuth theta
!      2) Rebuild global actuator points
!      3) Compute per-point forces (CPU or GPU)
!      4) MPI_Allreduce to accumulate tile contributions
!      5) Deposit smoothed forces on local F_turb
!
!  CALL:
!    call turbine_compute_forcing(F_turb, turbines, rho, u, v, w)
!--------------------------------------------------------------
subroutine turbine_forcing(F_turb_local, turbines_in, rho, u, v, w)
   use mod_dimensions, only : nx, ny, nz
   use mod_turbines
   use m_turbine_points
   use m_turbine_deposit
   use m_wtime
#ifdef MPI
   use m_mpi_decomp_init, only : j_start, j_end, mpi_rank
#endif
#ifdef _CUDA
   use cudafor
#endif
   implicit none

   ! Output forcing field on this tile
   real, intent(inout)            :: F_turb_local(3,0:nx+1,0:ny+1,0:nz+1)

   ! Turbine configuration (can be same as mod_turbines::turbines)
   type(turbine_t), intent(inout) :: turbines_in(:)

   ! Flow fields (tile-local storage; often device arrays under CUDA)
   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: w  (0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: rho,u,v,w
#endif

   integer :: np,i
   real, allocatable :: rho_h(:,:,:), u_h(:,:,:), v_h(:,:,:), w_h(:,:,:)

   call cpustart()

   ! 1. Update turbine azimuth
   turbines_in(:)%theta = turbines_in(:)%theta + turbines_in(:)%omegand

   ! 2. Construct global actuator points
   if (allocated(points_global)) deallocate(points_global)
   call turbine_distribute_points(turbines_in, points_global)
   np = size(points_global)

   if (np == 0) then
      F_turb_local = 0.0
      return
   end if

   ! 3. Allocate global force vectors in mod_turbines
   if (allocated(Fvec_local))  deallocate(Fvec_local)
   if (allocated(Fvec_global)) deallocate(Fvec_global)
   allocate(Fvec_local(3, np))
   allocate(Fvec_global(3, np))
   Fvec_local  = 0.0
   Fvec_global = 0.0
   call cpufinish(21)

   call cpustart()
#ifdef _CUDA
   ! 4. GPU path: compute point forces directly on device
   call turbine_point_forces_gpu(points_global, rho, u, v, w, Fvec_local, np)
#else
   ! 4. CPU path: copy fields and compute on host
   allocate(rho_h(0:nx+1,0:ny+1,0:nz+1))
   allocate(u_h  (0:nx+1,0:ny+1,0:nz+1))
   allocate(v_h  (0:nx+1,0:ny+1,0:nz+1))
   allocate(w_h  (0:nx+1,0:ny+1,0:nz+1))

   rho_h = rho
   u_h   = u
   v_h   = v
   w_h   = w

   call turbine_point_forces(points_global, rho_h, u_h, v_h, w_h, Fvec_local, np)

   deallocate(rho_h, u_h, v_h, w_h)
#endif
   call cpufinish(22)

   call cpustart()
#ifdef MPI
   ! 5. MPI reduction: sum contributions from all tiles
   call MPI_Allreduce(Fvec_local, Fvec_global, 3*np, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
#else
   Fvec_global = Fvec_local
#endif
   call cpufinish(23)

   call cpustart()
   ! 6. Clear local forcing field and deposit smoothed forces
   F_turb_local = 0.0

   call turbine_deposit(F_turb_local, points_global, Fvec_global, np)

   deallocate(Fvec_local, Fvec_global)
   deallocate(points_global)

   call cpufinish(24)
end subroutine turbine_forcing

end module m_turbine_forcing
