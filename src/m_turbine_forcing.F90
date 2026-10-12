

!=======================================================================
! Turbine forcing driver
!
! Coordinates the complete actuator-line / actuator-disk turbine
! forcing calculation.
!
! Computational sequence
! ----------------------
!
!   1. turbine_controller
!        Update turbine operating state (rotor speed, pitch and yaw).
!
!   2. turbine_distribute_points
!        Construct the actuator points from the current turbine state
!        and blade geometry.
!
!   3. turbine_point_kinematics_gpu
!        Interpolate flow variables to the actuator points and compute
!        axial velocity, relative tangential velocity and density.
!
!      For MPI runs the actuator-point kinematics are assembled on all
!      ranks using MPI_Allreduce before the blade forces are evaluated.
!
!   4. turbine_smearing_forces
!        Compute the aerodynamic actuator-point forces.  For ALM the
!        optional near-wake smearing correction modifies the velocity
!        used for the blade-element aerodynamic calculation.
!
!        Two independent switches select the aerodynamic corrections:
!
!           use_smearing   use_prandtl
!           ------------   -----------
!             .false.       .false.    basic ALM
!             .false.       .true.     ALM + Prandtl tip loss
!             .true.        .false.    ALM + smearing correction
!             .true.        .true.     ALM + smearing + Prandtl
!
!        The smearing correction is currently applied only to ALM.
!        Detailed aerodynamic and near-wake formulations are documented
!        in the corresponding force and smearing modules.
!
!   5. turbine_controller_feedback / turbine diagnostics
!        Calculate turbine loads, power and controller feedback from
!        the resulting actuator forces.
!
!   6. turbine_deposit_gpu
!        Spread the actuator-point forces onto the LBM grid using the
!        Gaussian force-deposition kernel and add them to
!        external_forcing.
!
! IMPORTANT:
!   The Gaussian smearing length used by the near-wake correction must
!   be identical to that used by the actuator-force deposition.
!
!=======================================================================
module m_turbine_forcing
contains
!      1) Compute local upstream wind and update turbine yaw
!      2) Update rotor azimuth theta
!      3) Rebuild global actuator point locations
!      4) Compute per-point forces (CPU or GPU)
!      5) MPI_Allreduce to accumulate tile contributions
!      6) Deposit smoothed forces on local external_forcing
subroutine turbine_forcing(external_forcing, turbines_in, rho, u, v, w, itimestep)
#ifdef MPI
   use mpi
#endif
   use mod_dimensions, only : nx, ny, nz
   use m_readinfile, only : p2l, dtcontrol, use_smearing, use_prandtl

   use mod_turbines, only : turbine_t, points_global
   use m_turbine_controller, only : turbine_controller, turbine_controller_feedback
   use m_turbine_distribute_points, only : turbine_distribute_points
   use m_turbine_point_kinematics_gpu, only : turbine_point_kinematics_gpu
   use m_turbine_smearing_forces, only : turbine_smearing_forces
   use m_turbine_diagnostics, only : turbine_diagnostics
   use m_turbine_diagnostics_timeseries, only : turbine_diagnostics_timeseries

   use m_turbine_deposit, only : turbine_deposit
   use m_turbine_deposit_gpu, only : turbine_deposit_gpu

   use m_wtime, only : cpustart, cpufinish
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
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
   real,            allocatable :: kin_local(:,:)    ! (3,np): uaxis, utheta, rho
   real,            allocatable :: kin_global(:,:)   ! global actuator-point kinematics
   real,            allocatable :: Fvec_global(:,:)  ! (3,np): corrected global forces

   real,            allocatable :: xg_h(:),yg_h(:),zg_h(:) ! (np)
#ifdef _CUDA
   real, device,    allocatable :: xg(:),yg(:),zg(:) ! (np)
#endif
#ifdef _CUDA
   real, device,    allocatable :: Fvec(:,:)         ! (3, np)
#else
   real,            allocatable :: Fvec(:,:)         ! (3, np)
#endif
#ifndef MPI
   integer :: mpi_rank=0
#endif
   real, parameter :: pi2=2.0*acos(-1.0)
   integer :: p, np, ierr, ncontrol, itimestep
   logical newfile

   call cpustart()
   ncontrol = max(1,nint(dtcontrol/p2l%time))
!-----------------------------------------------------------------------
! 1. Update turbine controllers:
!
!    - local upstream wind
!    - wind filtering
!    - yaw
!    - rotor RPM
!    - blade pitch
!-----------------------------------------------------------------------
   call turbine_controller(turbines_in,u,v,w,rho,itimestep)





!-----------------------------------------------------------------------
! 2. Update rotor azimuth.
!-----------------------------------------------------------------------
   turbines_in(:)%theta = modulo(turbines_in(:)%theta+turbines_in(:)%omegand,pi2)






! 3. Construct global actuator point locations and blade data stored in points_global(np)
   if (allocated(points_global)) deallocate(points_global)
   call turbine_distribute_points(turbines_in, points_global)
   np = size(points_global)

   call cpufinish(21)





   call cpustart()
! 4. Sample actuator-point kinematics on the owning MPI tile.
!
!    kin(1,p) = axial CFD velocity in rotor coordinates
!    kin(2,p) = relative tangential velocity omega*r-u.e_tan
!    kin(3,p) = interpolated density
!
!    The near-wake correction needs complete blades, so forces cannot be
!    finalized independently on each MPI tile at this stage.
   allocate(kin_local(3,np),kin_global(3,np),Fvec_global(3,np))
   kin_local  = 0.0
   kin_global = 0.0
   Fvec_global= 0.0

   call turbine_point_kinematics_gpu(points_global,rho,u,v,w,kin_local,np)
   call cpufinish(22)

#ifdef MPI
! 5. Assemble the complete actuator-point velocity field on every rank.
   call MPI_Allreduce(kin_local,kin_global,3*np,MPI_REAL,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
   kin_global=kin_local
#endif

! 6. Apply the near-wake Gaussian-smearing correction blade by blade and
!    compute the final actuator forces.  The correction is used for ALM
!    turbines only; ADM-R points use the uncorrected sampled velocities.
!
   call cpustart()
   call turbine_smearing_forces(turbines_in,points_global,kin_global, &
                                Fvec_global,np,p2l%length,p2l%vel, &
                                p2l%time,itimestep,use_smearing,use_prandtl)
   call cpufinish(23)



!-----------------------------------------------------------------------
! Update aerodynamic feedback for the next pitch-controller step.
!-----------------------------------------------------------------------
   call turbine_controller_feedback(turbines_in,points_global, &
                                    Fvec_global,np,itimestep)

! Turbine aerodynamic diagnostics
  if (mod(itimestep,10*ncontrol) == 0) then
      call turbine_diagnostics(turbines_in, points_global, Fvec_global, np)
  endif

!  if (mod(itimestep,1) == 0) then
     if (itimestep==1) newfile=.true.
     call turbine_diagnostics_timeseries(turbines_in,points_global,Fvec_global,np,itimestep,newfile)
     newfile=.false.
!  endif

! 6. Clear local forcing field and deposit smoothed forces
   call cpustart()
!-----------------------------------------------------------------------------------
#ifdef _CUDA
! device computation
      if (.not.allocated(xg_h)) allocate(xg_h(np))
      if (.not.allocated(yg_h)) allocate(yg_h(np))
      if (.not.allocated(zg_h)) allocate(zg_h(np))
      do p=1,np
         xg_h(p)=points_global(p)%xg
         yg_h(p)=points_global(p)%yg
         zg_h(p)=points_global(p)%zg
      enddo

      if (.not.allocated(xg)) allocate(xg(np))
      if (.not.allocated(yg)) allocate(yg(np))
      if (.not.allocated(zg)) allocate(zg(np))
      xg=xg_h
      yg=yg_h
      zg=zg_h
      if (.not.allocated(fvec)) allocate(fvec(3,np))
      fvec=fvec_global
      call turbine_deposit_gpu(external_forcing, xg, yg, zg, fvec, np)
      !print *,'sums =',mpi_rank,sum(external_forcing(:,1:nx,1:ny,1:nz)),sum(Fvec_global(:,:))
      !deallocate(fvec,xg,yg,zg)
!---------------------------------------------------------------------------------
#else
! host computation
      allocate(F_turb(3,0:nx+1,0:ny+1,0:nz+1))
      F_turb = 0.0
      call turbine_deposit(F_turb, points_global, Fvec_global, np)
      external_forcing=F_turb
      deallocate(F_turb)
#endif
!-----------------------------------------------------------------------------------

   deallocate(kin_local,kin_global,Fvec_global)
   deallocate(points_global)
   call cpufinish(24)

end subroutine turbine_forcing

end module m_turbine_forcing
