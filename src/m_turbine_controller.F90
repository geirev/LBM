module m_turbine_controller
   use mod_turbine_def, only : ratedrpm,ratedpower,pitch_min,pitch_max,pitchrate_max

   implicit none

!-----------------------------------------------------------------------
! Controller state.
!
! paero_last(n)        : most recent aerodynamic rotor power [W]
! torque_aero_last(n)  : most recent aerodynamic rotor torque [N m]
! pitch_integral(n)    : integral state of normalized rotor-speed error [s]
! power_initialized(n) : true after first aerodynamic torque measurement
!-----------------------------------------------------------------------
   real,    allocatable, save :: paero_last(:)
   real,    allocatable, save :: torque_aero_last(:)
   real,    allocatable, save :: pitch_integral(:)
   logical, allocatable, save :: power_initialized(:)

   logical, save :: controller_initialized = .false.

!-----------------------------------------------------------------------
! SWT-6.0-154 surrogate drivetrain/controller parameters.
!
! rotor_inertia is the combined rotor/direct-drive-generator inertia
! referred to the rotor [kg m^2].  A reliable public Siemens value is
! not available, so this is explicitly a tunable model parameter.
!
! rpm_min is the lower variable-speed operating limit.
!
! Region 2 generator torque is
!
!       Qgen = Kgen*omega^2
!
! where Kgen is calculated from ratedpower and ratedrpm so that Region 2
! joins continuously to rated power at rated rotor speed.
!-----------------------------------------------------------------------
   real, save :: rotor_inertia = 2.48e7      ! [kg m^2], initial estimate
   real, save :: rpm_min       = 4.0         ! [rev/min]

!-----------------------------------------------------------------------
! Above-rated pitch PI gains based on normalized rotor-speed error:
!
!       error = (omega-omega_rated)/omega_rated
!
! pitch_kp [deg]
! pitch_ki [deg/s]
!
! These are initial tuning values and may require adjustment after the
! drivetrain inertia and aerodynamic response have been calibrated.
!-----------------------------------------------------------------------
   real, save :: pitch_kp = 10.0
   real, save :: pitch_ki = 0.5

contains


!=======================================================================
! Turbine controller
!
! BELOW RATED (Region 2):
!
!    Qgen = Kgen*omega^2
!
!    J domega/dt = Qaero-Qgen
!
!    pitch -> pitch_min
!
! Rotor speed is therefore determined dynamically by the aerodynamic
! torque from the ALM and the generator torque.  No estimate of incoming
! wind speed is required for RPM control.
!
! ABOVE RATED (Region 3):
!
!    generator power is limited approximately to ratedpower
!
!    Qgen = ratedpower/omega
!
! and collective blade pitch regulates rotor speed toward ratedrpm.
!
! The rotor speed remains dynamic in Region 3; it is not hard-clamped to
! ratedrpm.  This allows gusts and turbulent/waked inflow to produce a
! physically meaningful drivetrain response.
!=======================================================================
subroutine turbine_controller(turbines_in,u,v,w,rho,itimestep)

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions, only : nx,ny,nz

   use mod_turbines, only : turbine_t, &
                            uavg_f,vavg_f,wavg_f, &
                            windfilter_initialized

   use m_readinfile, only : p2l,udir,localwind,dtcontrol,filter_time, &
                            yawrate_max,yaw_deadband

   use m_turbine_local_wind
   use m_turbine_yaw_controller

   implicit none

   type(turbine_t), intent(inout) :: turbines_in(:)

   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: w  (0:nx+1,0:ny+1,0:nz+1)

#ifdef _CUDA
   attributes(device) :: rho,u,v,w
#endif

   integer, intent(in) :: itimestep

   integer :: n,ncontrol

   real, parameter :: pi  = acos(-1.0)
   real, parameter :: pi2 = 2.0*pi

   real :: dtcontrol_actual
   real :: alpha

   real :: uavg,vavg,wavg
   real :: speed,winddir,unormal

   real :: omega
   real :: omega_min
   real :: omega_rated
   real :: kgen
   real :: torque_gen
   real :: domega


!-----------------------------------------------------------------------
! Initialize controller state.
!-----------------------------------------------------------------------
   call turbine_controller_initialize(size(turbines_in))


!-----------------------------------------------------------------------
! Controller interval.
!-----------------------------------------------------------------------
   ncontrol = max(1,nint(dtcontrol/p2l%time))

   dtcontrol_actual = real(ncontrol)*p2l%time

   alpha = dtcontrol_actual/(filter_time + dtcontrol_actual)


!-----------------------------------------------------------------------
! Controller acts only at the specified interval.
!-----------------------------------------------------------------------
   if (mod(itimestep,ncontrol) /= 0) return


!-----------------------------------------------------------------------
! Rotor-speed limits and Region-2 generator constant.
!-----------------------------------------------------------------------
   omega_min   = rpm_min*pi2/60.0
   omega_rated = ratedrpm*pi2/60.0

   if (omega_rated <= 0.0) then
      write(*,*) 'ERROR turbine_controller: ratedrpm must be positive'
      error stop
   endif

   if (ratedpower <= 0.0) then
      write(*,*) 'ERROR turbine_controller: ratedpower must be positive'
      error stop
   endif

   if (rotor_inertia <= 0.0) then
      write(*,*) 'ERROR turbine_controller: rotor_inertia must be positive'
      error stop
   endif

   kgen = ratedpower/omega_rated**3


!=======================================================================
! Turbine loop.
!=======================================================================
   do n = 1,size(turbines_in)


!-----------------------------------------------------------------------
! Local upstream wind is required only for local yaw control.
! RPM control no longer requires an estimate of incoming wind speed.
!-----------------------------------------------------------------------
      if (localwind == 1) then

         call turbine_local_wind(turbines_in(n),u,v,w,rho, &
                                 1.0,6, &
                                 uavg,vavg,wavg, &
                                 speed,winddir,unormal)


!-----------------------------------------------------------------------
! Low-pass filter local velocity components used by the yaw controller.
!-----------------------------------------------------------------------
         if (.not.windfilter_initialized(n)) then

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


!-----------------------------------------------------------------------
! Yaw controller.
!-----------------------------------------------------------------------
      if (localwind == 1 .and. windfilter_initialized(n)) then

         winddir = atan2(vavg_f(n),uavg_f(n))*360.0/pi2

         call turbine_yaw_controller(winddir,dtcontrol_actual,1, &
              turbines_in(n:n)%yaw,yawrate_max,yaw_deadband)

      elseif (localwind == 2) then

         turbines_in(n)%yaw = &
              wrap_180_controller(udir)*pi2/360.0

      endif


!-----------------------------------------------------------------------
! Current physical rotor angular velocity [rad/s].
!-----------------------------------------------------------------------
      omega = turbines_in(n)%omegand/p2l%time


!-----------------------------------------------------------------------
! Before the first aerodynamic feedback measurement, retain the existing
! rotor speed.  This avoids inventing an aerodynamic torque at startup.
!-----------------------------------------------------------------------
      if (.not.power_initialized(n)) then

         omega = max(omega_min,min(omega_rated,omega))

         call turbine_pitch_below_rated( &
              turbines_in(n)%pitchangle, &
              pitch_integral(n), &
              dtcontrol_actual)

      else


!-----------------------------------------------------------------------
! Generator torque.
!
! Below rated:
!
!       Qgen = Kgen*omega^2
!
! Above rated:
!
!       Qgen = Prated/omega
!
! The two expressions are identical at omega=omega_rated because
!
!       Kgen = Prated/omega_rated^3.
!-----------------------------------------------------------------------
         if (omega < omega_rated) then

            torque_gen = kgen*omega**2

         else

            torque_gen = ratedpower/max(omega,1.0e-12)

         endif


!-----------------------------------------------------------------------
! Rotor/direct-drive-generator dynamics.
!
!       J domega/dt = Qaero-Qgen
!-----------------------------------------------------------------------
         domega = (torque_aero_last(n)-torque_gen)/rotor_inertia

         omega = omega + domega*dtcontrol_actual


!-----------------------------------------------------------------------
! Lower rotor-speed limit.
!
! The upper speed is deliberately not hard-clamped to ratedrpm.  Above
! rated, the pitch controller acts on rotor-speed error and brings the
! rotor back toward rated speed dynamically.
!-----------------------------------------------------------------------
         omega = max(omega_min,omega)


!-----------------------------------------------------------------------
! Pitch control.
!
! Below rated: return toward pitch_min.
! Above rated: regulate rotor speed toward ratedrpm.
!-----------------------------------------------------------------------
         if (omega < omega_rated) then

            call turbine_pitch_below_rated( &
                 turbines_in(n)%pitchangle, &
                 pitch_integral(n), &
                 dtcontrol_actual)

         else

            call turbine_pitch_speed_controller( &
                 omega,omega_rated, &
                 dtcontrol_actual, &
                 turbines_in(n)%pitchangle, &
                 pitch_integral(n))

         endif

      endif


!-----------------------------------------------------------------------
! Store rotor angular increment in lattice units [rad/timestep].
!-----------------------------------------------------------------------
      turbines_in(n)%omegand = omega*p2l%time

   enddo

end subroutine turbine_controller


!=======================================================================
! Above-rated collective pitch controller.
!
! Uses normalized rotor-speed error:
!
!       error = (omega-omega_rated)/omega_rated
!
! Thus:
!
!       error > 0   -> increase pitch
!       error < 0   -> decrease pitch
!
! The PI controller is
!
!       beta_target = beta_min + Kp*error + Ki*integral(error dt)
!
! A pitch-rate limit and simple integral anti-windup are applied.
!=======================================================================
subroutine turbine_pitch_speed_controller(omega,omega_rated,dt,pitch,integral_error)

   implicit none

   real, intent(in)    :: omega
   real, intent(in)    :: omega_rated
   real, intent(in)    :: dt

   real, intent(inout) :: pitch
   real, intent(inout) :: integral_error

   real :: error
   real :: integral_trial

   real :: pitch_target
   real :: pitch_trial

   real :: dpitch
   real :: dpitch_max


   if (omega_rated <= 0.0) return


!-----------------------------------------------------------------------
! Normalized rotor-speed error [-].
!-----------------------------------------------------------------------
   error = (omega-omega_rated)/omega_rated


!-----------------------------------------------------------------------
! Trial integral.
!-----------------------------------------------------------------------
   integral_trial = integral_error + error*dt


!-----------------------------------------------------------------------
! Trial PI output [deg].
!-----------------------------------------------------------------------
   pitch_trial = pitch_min + &
                 pitch_kp*error + &
                 pitch_ki*integral_trial


!-----------------------------------------------------------------------
! Integral anti-windup.
!
! Only accept integration while the unconstrained controller output
! remains inside the allowed pitch range.
!-----------------------------------------------------------------------
   if (pitch_trial >= pitch_min .and. &
       pitch_trial <= pitch_max) then

      integral_error = integral_trial

   endif


!-----------------------------------------------------------------------
! Final target pitch.
!-----------------------------------------------------------------------
   pitch_target = pitch_min + &
                  pitch_kp*error + &
                  pitch_ki*integral_error

   pitch_target = min(pitch_max,max(pitch_min,pitch_target))


!-----------------------------------------------------------------------
! Pitch-rate limiter [deg].
!-----------------------------------------------------------------------
   dpitch_max = pitchrate_max*dt

   dpitch = pitch_target-pitch

   if (abs(dpitch) > dpitch_max) then
      dpitch = sign(dpitch_max,dpitch)
   endif

   pitch = pitch + dpitch

   pitch = min(pitch_max,max(pitch_min,pitch))

end subroutine turbine_pitch_speed_controller


!=======================================================================
! Below-rated pitch controller.
!
! Below rated, the blades are returned toward pitch_min while rotor
! speed is determined by aerodynamic and generator torque.
!=======================================================================
subroutine turbine_pitch_below_rated(pitch,integral_error,dt)

   implicit none

   real, intent(inout) :: pitch
   real, intent(inout) :: integral_error

   real, intent(in) :: dt

   real :: dpitch
   real :: dpitch_max


!-----------------------------------------------------------------------
! Reset integral state when leaving the above-rated region.
!-----------------------------------------------------------------------
   integral_error = 0.0


!-----------------------------------------------------------------------
! Return pitch gradually toward pitch_min.
!-----------------------------------------------------------------------
   dpitch = pitch_min-pitch

   dpitch_max = pitchrate_max*dt

   if (abs(dpitch) > dpitch_max) then
      dpitch = sign(dpitch_max,dpitch)
   endif

   pitch = pitch + dpitch

   pitch = min(pitch_max,max(pitch_min,pitch))

end subroutine turbine_pitch_below_rated


!=======================================================================
! Aerodynamic feedback update.
!
! Must be called after Fvec_global has been formed by MPI_Allreduce.
!
! For each turbine:
!
!       Q = sum[(r x F) . e_axis]       [LB torque]
!
!       Qphys = Q * torque_conv         [N m]
!
!       Paero = Qphys * omega           [W]
!
! Both aerodynamic torque and power are stored for the next controller
! update.  The torque is used directly in the drivetrain equation.
!=======================================================================
subroutine turbine_controller_feedback(turbines_in,points_global, &
                                       Fvec_global,np,itimestep)

   use mod_turbines, only : turbine_t,point_t
   use m_readinfile, only : p2l,dtcontrol
   use m_turbine_rotor_basis

   implicit none

   type(turbine_t), intent(in) :: turbines_in(:)
   type(point_t),   intent(in) :: points_global(:)

   integer, intent(in) :: np
   integer, intent(in) :: itimestep

   real, intent(in) :: Fvec_global(3,np)

   integer :: n,p
   integer :: nturb
   integer :: ncontrol

   real :: e_axis(3),e1(3),e2(3)

   real :: rvec(3)
   real :: fpoint(3)
   real :: torque_vec(3)

   real :: tq
   real :: torque

   real :: torque_conv
   real :: torque_phys
   real :: omega


!-----------------------------------------------------------------------
! Use the actual turbine-array size rather than the external nturbines
! variable. This guarantees consistency with turbines_in(:).
!-----------------------------------------------------------------------
   nturb = size(turbines_in)

   call turbine_controller_initialize(nturb)


!-----------------------------------------------------------------------
! Basic consistency checks.
!-----------------------------------------------------------------------
   if (np /= size(points_global)) then
      write(*,*) 'ERROR turbine_controller_feedback:'
      write(*,*) 'np                  = ',np
      write(*,*) 'size(points_global) = ',size(points_global)
      error stop
   endif

   if (size(Fvec_global,1) /= 3) then
      write(*,*) 'ERROR turbine_controller_feedback:'
      write(*,*) 'size(Fvec_global,1) = ',size(Fvec_global,1)
      error stop
   endif

   if (size(Fvec_global,2) /= np) then
      write(*,*) 'ERROR turbine_controller_feedback:'
      write(*,*) 'np                  = ',np
      write(*,*) 'size(Fvec_global,2)  = ',size(Fvec_global,2)
      error stop
   endif

   if (size(paero_last) /= nturb .or. &
       size(torque_aero_last) /= nturb .or. &
       size(power_initialized) /= nturb) then

      write(*,*) 'ERROR turbine_controller_feedback: controller arrays'
      write(*,*) 'nturb                       = ',nturb
      write(*,*) 'size(paero_last)             = ',size(paero_last)
      write(*,*) 'size(torque_aero_last)       = ',size(torque_aero_last)
      write(*,*) 'size(power_initialized)      = ',size(power_initialized)
      error stop

   endif


!-----------------------------------------------------------------------
! Only update aerodynamic feedback at controller intervals.
!-----------------------------------------------------------------------
   ncontrol = max(1,nint(dtcontrol/p2l%time))

   if (mod(itimestep,ncontrol) /= 0) return


!-----------------------------------------------------------------------
! Lattice torque -> physical torque:
!
!       Qphys = QLB * rho_c * L_c^5 / T_c^2
!
! [N m]
!-----------------------------------------------------------------------
   torque_conv = p2l%rho*p2l%length**5/p2l%time**2


!=======================================================================
! Turbine loop.
!=======================================================================
   do n = 1,nturb

      torque = 0.0


!-----------------------------------------------------------------------
! Rotor basis.
!-----------------------------------------------------------------------
      call turbine_rotor_basis(turbines_in(n)%yaw, &
                               turbines_in(n)%tilt, &
                               e_axis,e1,e2)


!-----------------------------------------------------------------------
! Integrate aerodynamic torque.
!-----------------------------------------------------------------------
      do p = 1,np

         if (points_global(p)%iturb /= n) cycle

         fpoint(:) = Fvec_global(:,p)

         rvec(1) = points_global(p)%xg - turbines_in(n)%xhub
         rvec(2) = points_global(p)%yg - turbines_in(n)%yhub
         rvec(3) = points_global(p)%zg - turbines_in(n)%zhub


!-----------------------------------------------------------------------
! r x F
!-----------------------------------------------------------------------
         torque_vec(1) = rvec(2)*fpoint(3) - &
                         rvec(3)*fpoint(2)

         torque_vec(2) = rvec(3)*fpoint(1) - &
                         rvec(1)*fpoint(3)

         torque_vec(3) = rvec(1)*fpoint(2) - &
                         rvec(2)*fpoint(1)


!-----------------------------------------------------------------------
! Q = (r x F) . e_axis
!-----------------------------------------------------------------------
         tq = dot_product(torque_vec,e_axis)

         torque = torque+tq

      enddo


!-----------------------------------------------------------------------
! Physical aerodynamic torque [N m].
!-----------------------------------------------------------------------
      torque_phys = torque*torque_conv

      torque_aero_last(n) = torque_phys


!-----------------------------------------------------------------------
! Physical angular velocity [rad/s].
!-----------------------------------------------------------------------
      omega = turbines_in(n)%omegand/p2l%time


!-----------------------------------------------------------------------
! Aerodynamic mechanical power [W].
!-----------------------------------------------------------------------
      paero_last(n) = torque_phys*omega

      power_initialized(n) = .true.

   enddo

end subroutine turbine_controller_feedback


!=======================================================================
! Allocate or resize controller state.
!=======================================================================
subroutine turbine_controller_initialize(nturb)

   implicit none

   integer, intent(in) :: nturb


   if (nturb <= 0) then
      write(*,*) 'ERROR turbine_controller_initialize: nturb = ',nturb
      error stop
   endif


!-----------------------------------------------------------------------
! If already allocated with the correct size, nothing needs to be done.
!-----------------------------------------------------------------------
   if (controller_initialized) then

      if (allocated(paero_last) .and. &
          allocated(torque_aero_last) .and. &
          allocated(pitch_integral) .and. &
          allocated(power_initialized)) then

         if (size(paero_last)         == nturb .and. &
             size(torque_aero_last)   == nturb .and. &
             size(pitch_integral)     == nturb .and. &
             size(power_initialized)  == nturb) return

      endif

   endif


!-----------------------------------------------------------------------
! Allocate/reallocate controller state.
!-----------------------------------------------------------------------
   if (allocated(paero_last))         deallocate(paero_last)
   if (allocated(torque_aero_last))   deallocate(torque_aero_last)
   if (allocated(pitch_integral))     deallocate(pitch_integral)
   if (allocated(power_initialized))  deallocate(power_initialized)

   allocate(paero_last(nturb))
   allocate(torque_aero_last(nturb))
   allocate(pitch_integral(nturb))
   allocate(power_initialized(nturb))

   paero_last         = 0.0
   torque_aero_last   = 0.0
   pitch_integral     = 0.0
   power_initialized  = .false.

   controller_initialized = .true.

end subroutine turbine_controller_initialize


!=======================================================================
! Wrap angle to [-180,180) degrees.
!=======================================================================
pure real function wrap_180_controller(angle)

   implicit none

   real, intent(in) :: angle

   wrap_180_controller = modulo(angle+180.0,360.0)-180.0

end function wrap_180_controller


end module m_turbine_controller
