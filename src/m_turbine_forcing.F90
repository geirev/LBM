module m_turbine_forcing
contains
!      1) Update rotor azimuth theta
!      2) Rebuild global actuator point locations
!      3) Compute per-point forces (CPU or GPU)
!      4) MPI_Allreduce to accumulate tile contributions
!      5) Deposit smoothed forces on local F_turb
subroutine turbine_forcing(external_forcing, turbines_in, rho, u, v, w)
#ifdef MPI
   use mpi
#endif
   use mod_dimensions, only : nx, ny, nz
   use mod_turbines, only : turbine_t,point_t,points_global

   use m_turbine_distribute_points
   use m_turbine_point_forces_gpu
   use m_turbine_point_forces

   use m_turbine_deposit

   use m_turbines_bounding_box
   use m_wtime
#ifdef MPI
   use m_mpi_decomp_init, only : j_start, j_end, mpi_rank
#endif
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   real,            intent(inout) :: external_forcing(3,0:nx+1,0:ny+1,0:nz+1) ! Output forcing field on this tile
   type(turbine_t), intent(inout) :: turbines_in(:)                 ! Turbine configuration
   real,            intent(in)    :: rho(0:nx+1,0:ny+1,0:nz+1)
   real,            intent(in)    :: u  (0:nx+1,0:ny+1,0:nz+1)
   real,            intent(in)    :: v  (0:nx+1,0:ny+1,0:nz+1)
   real,            intent(in)    :: w  (0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: rho,u,v,w,external_forcing
#endif

   real,            allocatable :: F_turb(:,:,:,:)
   real,            allocatable :: Fvec_local(:,:)   ! (3, np)
   real,            allocatable :: Fvec_global(:,:)  ! (3, np)


   integer :: np,i,ierr
   real, allocatable :: rho_h(:,:,:), u_h(:,:,:), v_h(:,:,:), w_h(:,:,:)
   integer krad
#ifndef MPI
   integer :: mpi_rank=0
#endif

   call cpustart()
! 1. Update turbine azimuth
!  if (tipspeed /= 0.0) then
!     compute new turbrpm for each turbine
!     do n=1,nrturbines
!        turbrpm=f(u,tipspeedratio)
!        omega = pi2 * turbrpm / 60.0
!        turbines(n)%omegand    = omega * p2l%time
!     enddo
!  endif

   turbines_in(:)%theta = turbines_in(:)%theta + turbines_in(:)%omegand

! 2. Construct global actuator point locations and blade data stored in points_global(np)
   if (allocated(points_global)) deallocate(points_global)
   call turbine_distribute_points(turbines_in, points_global)
   np = size(points_global)
   call cpufinish(21)

   call cpustart()
! 3. Allocate global force vectors
   if (allocated(Fvec_local))  deallocate(Fvec_local) ; allocate(Fvec_local(3, np))
   if (allocated(Fvec_global)) deallocate(Fvec_global); allocate(Fvec_global(3, np))
   Fvec_local  = 0.0
   Fvec_global = 0.0

! 4. Compute point forces directly on device
   call turbine_point_forces_gpu(points_global, rho, u, v, w, Fvec_local, np)
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
   allocate(F_turb(3,0:nx+1,0:ny+1,0:nz+1))
   F_turb = 0.0
   call turbine_deposit(F_turb, points_global, Fvec_global, np, krad)
!   call turbines_bounding_box(points_global, np, krad)
   external_forcing=F_turb

   deallocate(Fvec_local, Fvec_global)
   deallocate(points_global)
!   deallocate(F_turb)

   call cpufinish(24)
end subroutine turbine_forcing

end module m_turbine_forcing
