!-----------------------------------------------------------------------
! Turbine forcing and yaw control
!-----------------------------------------------------------------------
! This routine computes the actuator-line turbine forcing and updates the
! yaw angle of each turbine using the locally resolved upstream wind.
!
! The main sequence is:
!
!   1. Local wind measurement, yaw control and rotor-speed control
!   2. Update rotor azimuth
!   3. Rebuild global actuator-point locations
!   4. Compute actuator-point forces
!   5. Combine point forces across MPI tiles
!   6. Deposit the smoothed turbine forces onto the local CFD grid!
!
!-----------------------------------------------------------------------
! Local wind, yaw control and tip-speed-ratio control
!-----------------------------------------------------------------------
!
! The local upstream wind is sampled every dtcontrol seconds using
! turbine_local_wind(). The same filtered velocity components are used
! for both yaw control and, when requested, rotor-speed control.
!
! A local wind measurement is performed when either
!
!      localwind == 1
!
! or
!
!      tipspeedratio > 0
!
! so TSR control can be used independently of the local yaw controller.
!
! The measured velocity components are low-pass filtered according to
!
!      q_f = q_f + alpha*(q-q_f)
!
! where
!
!      alpha = dtcontrol_actual/(filter_time + dtcontrol_actual).
!
! On the first measurement the filtered velocity is initialized directly
! from the measured velocity to avoid a startup transient.
!
! Yaw control:
!
!   localwind = 0 : use yaw values specified in infile.in
!   localwind = 1 : yaw each turbine toward its filtered local wind
!   localwind = 2 : set yaw from the externally imposed direction udir
!
! Rotor-speed control:
!
!   tipspeedratio <= 0 :
!      use the fixed RPM specified in infile.in.
!
!   tipspeedratio > 0 :
!      override the fixed RPM individually for each turbine using the
!      magnitude of its filtered local upstream velocity.
!
! The target tip-speed ratio is
!
!      lambda = omega*R/U
!
! giving
!
!      omega = lambda*U/R
!
! where R = hubradius + rotorradius is the physical rotor radius [m]
! and U is the filtered local upstream wind speed [m/s].
!
! The resulting angular velocity is stored in nondimensional form as
!
!      omegand = omega*p2l%time
!
! and is subsequently used to advance the rotor azimuth:
!
!      theta = theta + omegand
!
! Thus, when tipspeedratio > 0, turbines operating in different local
! wind conditions can rotate at different RPM while maintaining the
! same prescribed target tip-speed ratio.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 1. Local wind measurement and yaw control
!-----------------------------------------------------------------------
! The yaw controller uses the locally resolved wind direction rather than
! the externally imposed inflow direction.
!
! The controller is executed approximately every dtcontrol seconds:
!
!      ncontrol        = nint(dtcontrol/p2l%time)
!      dtcontrol_actual = ncontrol*p2l%time
!
! The actual controller interval dtcontrol_actual is used in both the
! low-pass filter and the yaw-rate limiter.
!
! For each turbine, turbine_local_wind() evaluates the velocity over a
! circular sampling disk:
!
!   - centered on the turbine axis,
!   - parallel to the current rotor plane,
!   - located 1 rotor diameter upstream,
!   - with diameter equal to the rotor diameter,
!   - discretized using nsrad=6 radial sampling intervals.
!
! The routine returns the disk-averaged velocity components
!
!      uavg, vavg, wavg
!
! together with the mean speed, instantaneous horizontal wind direction,
! and velocity normal to the rotor plane.
!
! The local measurement therefore accounts for spatial variations in the
! resolved flow caused by wakes, blockage, and other turbine interactions.
!
!-----------------------------------------------------------------------
! Low-pass filtering
!-----------------------------------------------------------------------
! The disk-averaged velocity components are filtered independently using
! a first-order exponential low-pass filter:
!
!      q_f = q_f + alpha*(q-q_f)
!
! with
!
!      alpha = dtcontrol_actual/(filter_time + dtcontrol_actual)
!
! where filter_time is the filter time constant in physical seconds.
!
! The current setting is
!
!      filter_time = 5.0 s
!
! On the first controller update, the filtered velocity is initialized
! directly from the first local-wind measurement to avoid a transient
! associated with initialization from zero.
!
! The wind direction is calculated AFTER filtering:
!
!      winddir = atan2(vavg_f,uavg_f)
!
! This avoids angular averaging and the associated discontinuity at
! +/-180 or 0/360 degrees.
!
! winddir is converted from radians to degrees before being passed to
! turbine_yaw_controller().
!
!-----------------------------------------------------------------------
! Yaw controller
!-----------------------------------------------------------------------
! turbine_yaw_controller() changes each turbine yaw angle subject to:
!
!      yawrate_max  = 0.3 deg/s
!      yaw_deadband = 5.0 deg
!
! The controller therefore:
!
!   - ignores yaw errors within the deadband,
!   - limits yaw motion to yawrate_max,
!   - aligns each turbine gradually with its filtered local wind direction.
!
! Internally, turbine yaw is stored in radians, while winddir,
! yawrate_max, and yaw_deadband are supplied to the controller in degrees.
!
!-----------------------------------------------------------------------
! MPI/GPU treatment of local wind
!-----------------------------------------------------------------------
! The macroscopic fields rho,u,v,w remain on the GPU.
!
! Each MPI tile evaluates only the part of the upstream sampling disk
! belonging to its local j-domain. The local velocity sums and number of
! valid sampling points are combined with MPI_Allreduce so that all MPI
! ranks obtain the same disk-averaged wind for each turbine.
!
! Consequently, each MPI rank obtains the same yaw update and the turbine
! geometry remains synchronized across tiles.
!
!-----------------------------------------------------------------------
! 2. Rotor azimuth
!-----------------------------------------------------------------------
! The blade azimuth is advanced according to the nondimensional angular
! velocity stored for each turbine:
!
!      theta = theta + omegand
!
!-----------------------------------------------------------------------
! 3. Actuator-point locations
!-----------------------------------------------------------------------
! turbine_distribute_points() rebuilds the actuator-point coordinates
! using the current turbine yaw, tilt, azimuth, and blade geometry.
!
!-----------------------------------------------------------------------
! 4. Actuator-point forces
!-----------------------------------------------------------------------
! turbine_point_forces_gpu() interpolates the local flow at each actuator
! point and computes the corresponding blade-force vector.
!
! Under MPI, each tile contributes only for actuator points whose
! interpolation stencil belongs to that tile.
!
!-----------------------------------------------------------------------
! 5. MPI force reduction
!-----------------------------------------------------------------------
! MPI_Allreduce sums the partial actuator-point forces from all MPI tiles,
! producing one global force vector for every actuator point.
!
!-----------------------------------------------------------------------
! 6. Force deposition
!-----------------------------------------------------------------------
! The global actuator-point forces are distributed onto the CFD grid using
! the turbine-force deposition routine. Under CUDA this operation is
! performed directly on the GPU.
!
!-----------------------------------------------------------------------
! Alternative yaw specification
!-----------------------------------------------------------------------
! For testing, local yaw control can be bypassed and the turbine yaw can
! instead be aligned directly with the externally imposed wind direction:
!
!      turbines_in(n)%yaw = (wrap_180(udir)/360.0)*pi2
!
! This is useful for comparison with the local-wind yaw controller.
!-----------------------------------------------------------------------

module m_turbine_forcing
contains
pure real function wrap_180(angle)
   implicit none
   real, intent(in) :: angle

   wrap_180 = modulo(angle + 180.0, 360.0) - 180.0

end function wrap_180

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
   use mod_dimensions, only : nx, ny, nz, nyg
   use m_readinfile, only : uini,udir,nturbines,p2l,localwind,dtcontrol,yawrate_max,yaw_deadband,filter_time,tipspeedratio

   use mod_turbines, only : turbine_t,point_t,points_global, uavg_f, vavg_f, wavg_f,windfilter_initialized
   use mod_turbine_def, only : hubradius,rotorradius,nrchords,relm

   use m_turbine_local_wind
   use m_turbine_yaw_controller
   use m_turbine_distribute_points
   use m_turbine_point_forces_gpu
   use m_turbine_diagnostics

   use m_turbine_deposit
   use m_turbine_deposit_gpu

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

   real,            allocatable :: xg_h(:),yg_h(:),zg_h(:) ! (np)
#ifdef _CUDA
   real, device,    allocatable :: xg(:),yg(:),zg(:) ! (np)
#endif
#ifdef _CUDA
   real, device,    allocatable :: Fvec(:,:)         ! (3, np)
#else
   real,            allocatable :: Fvec(:,:)         ! (3, np)
#endif
   logical :: lgpu=.true.

   integer :: p,np,i,ierr
   real, allocatable :: rho_h(:,:,:), u_h(:,:,:), v_h(:,:,:), w_h(:,:,:)
#ifndef MPI
   integer :: mpi_rank=0
#endif
   real, parameter :: pi=acos(-1.0)
   real, parameter :: pi2=2.0*pi
   integer n,ncontrol,itimestep
   real uavg, vavg, wavg, speed, winddir, unormal
   real dtcontrol_actual,alpha
   real :: radius, ulocal, omega, rpm


   call cpustart()


! Target controller interval in seconds and number of time steps per approximate second
   ncontrol = max(1,nint(dtcontrol/p2l%time))

! Actual controller interval in seconds represented by ncontrol model steps
   dtcontrol_actual = real(ncontrol)*p2l%time


   alpha = dtcontrol_actual/(filter_time + dtcontrol_actual)

   if (mod(itimestep,ncontrol) == 0) then

      do n = 1, nturbines

         !---------------------------------------------------------
         ! Local upstream wind is needed for:
         !   localwind = 1       : local yaw controller
         !   tipspeedratio > 0   : local rotor-speed controller
         !---------------------------------------------------------
         if (localwind == 1 .or. tipspeedratio > 0.0) then

            call turbine_local_wind(turbines_in(n),u,v,w,rho, &
                                    1.0,6, &
                                    uavg,vavg,wavg,speed,winddir,unormal)

            if (.not. windfilter_initialized(n)) then

               uavg_f(n) = uavg
               vavg_f(n) = vavg
               wavg_f(n) = wavg

               windfilter_initialized(n) = .true.

            else

               uavg_f(n) = uavg_f(n) + alpha*(uavg-uavg_f(n))
               vavg_f(n) = vavg_f(n) + alpha*(vavg-vavg_f(n))
               wavg_f(n) = wavg_f(n) + alpha*(wavg-wavg_f(n))

            endif

         endif


         !---------------------------------------------------------
         ! Yaw control
         !---------------------------------------------------------
         if (localwind == 1) then

            winddir = atan2(vavg_f(n),uavg_f(n))*360.0/pi2

            call turbine_yaw_controller(winddir,dtcontrol_actual,1, &
                 turbines_in(n)%yaw,yawrate_max,yaw_deadband)

         elseif (localwind == 2) then

            turbines_in(n)%yaw = (wrap_180(udir)/360.0)*pi2

         endif


         !---------------------------------------------------------
         ! Rotor speed from prescribed tip-speed ratio
         !
         ! lambda = omega R / U
         !---------------------------------------------------------
         if (tipspeedratio > 0.0 .and. windfilter_initialized(n)) then

            radius = rotorradius + hubradius       ! [m]

            ! Filtered local incoming wind speed [m/s]
            ulocal = sqrt(uavg_f(n)**2 + &
                          vavg_f(n)**2 + &
                          wavg_f(n)**2) * p2l%vel

            omega = tipspeedratio * ulocal / radius

            turbines_in(n)%omegand = omega * p2l%time

            ! Diagnostic physical RPM
            rpm = omega * 60.0/pi2

         endif

      enddo

   endif


   turbines_in(:)%theta = modulo(turbines_in(:)%theta + turbines_in(:)%omegand, pi2)

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

! Turbine aerodynamic diagnostics
   if (mod(itimestep,10*ncontrol) == 0) then
      call turbine_diagnostics(turbines_in, points_global, Fvec_global, np)
   endif

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

   deallocate(Fvec_local, Fvec_global)
   deallocate(points_global)
   call cpufinish(24)

end subroutine turbine_forcing

end module m_turbine_forcing
