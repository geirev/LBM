module m_turbine_controller
   use mod_turbine_def, only : hubradius,rotorradius, ratedrpm, ratedpower, pitch_min, pitch_max, pitchrate_max

   implicit none

!-----------------------------------------------------------------------
! Controller state.
!
! paero_last(n)      : most recently computed aerodynamic rotor power [W]
! power_initialized  : true after first aerodynamic power measurement
! pitch_integral     : integral state of normalized power error [s]
!-----------------------------------------------------------------------
   real,    allocatable, save :: paero_last(:)
   real,    allocatable, save :: pitch_integral(:)
   logical, allocatable, save :: power_initialized(:)

   logical, save :: controller_initialized = .false.

! PI gains based on normalized power error:
!
!    error = (Paero-Prated)/Prated
!
! Therefore:
!
!    pitch_kp [deg]
!    pitch_ki [deg/s]
!
! These are initial tuning values only.
   real, save :: pitch_kp = 10.0
   real, save :: pitch_ki = 0.5

contains


!=======================================================================
! Turbine controller
!
! Operating regions
! -----------------
!
! BELOW RATED:
!
!    omega = lambda_target*U/R
!    pitch -> pitch_min
!
! Thus the turbine maintains the prescribed optimal TSR.
!
! ABOVE RATED:
!
!    omega = omega_rated
!
! and collective blade pitch is controlled using the aerodynamic
! power measured on the previous controller update:
!
!    e = (Paero-Prated)/Prated
!
!    beta_target = beta_min + Kp*e + Ki*integral(e dt)
!
! Increasing pitch reduces angle of attack through
!
!    alpha = phi - twist - pitch
!
! which modifies Cl, Cd, torque and aerodynamic power.
!
! NOTE:
! The present turbine model prescribes rotor speed directly; it does
! not contain rotor/generator inertia. Therefore pitch regulates
! aerodynamic power above rated while RPM is explicitly capped.
!=======================================================================
subroutine turbine_controller(turbines_in,u,v,w,rho,itimestep)

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions, only : nx,ny,nz

   use mod_turbines, only : turbine_t, &
                            uavg_f,vavg_f,wavg_f, &
                            windfilter_initialized


   use m_readinfile, only : p2l,udir,nturbines, &
                            localwind,dtcontrol,filter_time, &
                            yawrate_max,yaw_deadband,tipspeedratio

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

   real :: radius
   real :: ulocal

   real :: omega_tsr
   real :: omega_rated
   real :: omega

   real :: rpm_tsr


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


!=======================================================================
! Turbine loop.
!=======================================================================
   do n = 1,size(turbines_in)


!-----------------------------------------------------------------------
! Local upstream wind.
!
! Required for:
!
!    localwind == 1       : local yaw control
!    tipspeedratio > 0    : TSR/RPM control
!-----------------------------------------------------------------------
      if (localwind == 1 .or. tipspeedratio > 0.0) then

         call turbine_local_wind(turbines_in(n),u,v,w,rho, &
                                 1.0,6, &
                                 uavg,vavg,wavg, &
                                 speed,winddir,unormal)


!-----------------------------------------------------------------------
! Low-pass filter local velocity components.
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
! Rotor-speed and pitch control.
!-----------------------------------------------------------------------
      if (tipspeedratio > 0.0 .and. windfilter_initialized(n)) then

         radius = rotorradius + hubradius                   ! [m]

         ! Filtered local incoming wind speed [m/s].
         ulocal = sqrt(uavg_f(n)**2 + &
                       vavg_f(n)**2 + &
                       wavg_f(n)**2)*p2l%vel

         if (ulocal > 1.0e-12) then

!-----------------------------------------------------------------------
! Rotor speed required to maintain target TSR:
!
!       omega_TSR = lambda U/R
!-----------------------------------------------------------------------
            omega_tsr = tipspeedratio*ulocal/radius         ! [rad/s]

            rpm_tsr = omega_tsr*60.0/pi2                    ! [rev/min]


!-----------------------------------------------------------------------
! Rated rotor speed.
!-----------------------------------------------------------------------
            omega_rated = ratedrpm*pi2/60.0                 ! [rad/s]


!-----------------------------------------------------------------------
! Below-rated operation:
!
! Maintain the prescribed optimal TSR until either the TSR-based rotor
! speed reaches rated RPM AND the aerodynamic power is approaching
! rated power.
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! Rotor-speed and pitch control.
!-----------------------------------------------------------------------
            if (.not.power_initialized(n)) then

               ! No aerodynamic feedback available yet.
               omega = min(omega_tsr,omega_rated)

               call turbine_pitch_below_rated( &
                    turbines_in(n)%pitchangle, &
                    pitch_integral(n), &
                    dtcontrol_actual)

            elseif (rpm_tsr <= ratedrpm) then

               !------------------------------------------------------------
               ! Region 2:
               ! Maintain prescribed optimal TSR.
               !------------------------------------------------------------
               omega = omega_tsr

               call turbine_pitch_below_rated( &
                    turbines_in(n)%pitchangle, &
                    pitch_integral(n), &
                    dtcontrol_actual)

            else

               !------------------------------------------------------------
               ! Region 3:
               ! Rotor has reached rated RPM.
               !
               ! Hold RPM at rated value and regulate aerodynamic power
               ! using collective blade pitch.
               !------------------------------------------------------------
               omega = omega_rated

               call turbine_pitch_power_controller( &
                    paero_last(n),ratedpower, &
                    dtcontrol_actual, &
                    turbines_in(n)%pitchangle, &
                    pitch_integral(n))

            endif
!-----------------------------------------------------------------------
! Store rotor angular increment in lattice units [rad/timestep].
!-----------------------------------------------------------------------
            turbines_in(n)%omegand = omega*p2l%time

         endif

      endif

   enddo

end subroutine turbine_controller


!=======================================================================
! Above-rated collective pitch controller.
!
! Uses normalized aerodynamic power error:
!
!      error = (Paero-Prated)/Prated
!
! Thus:
!
!      error > 0   -> increase pitch
!      error < 0   -> decrease pitch
!
! The PI controller is
!
!      beta_target =
!           beta_min + Kp*error + Ki*integral(error dt)
!
! A pitch-rate limit and simple integral anti-windup are applied.
!=======================================================================
subroutine turbine_pitch_power_controller(paero,prated,dt,pitch,integral_error)

   implicit none

   real, intent(in)    :: paero              ! [W]
   real, intent(in)    :: prated             ! [W]
   real, intent(in)    :: dt                 ! [s]

   real, intent(inout) :: pitch              ! [deg]
   real, intent(inout) :: integral_error     ! [s]

   real :: error
   real :: integral_trial

   real :: pitch_target
   real :: pitch_trial

   real :: dpitch
   real :: dpitch_max


   if (prated <= 0.0) return


!-----------------------------------------------------------------------
! Normalized power error [-].
!-----------------------------------------------------------------------
   error = (paero-prated)/prated


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

end subroutine turbine_pitch_power_controller


!=======================================================================
! Below-rated pitch controller.
!
! Below rated, the blades are returned toward pitch_min while rotor
! speed is controlled to maintain the target TSR.
!=======================================================================
subroutine turbine_pitch_below_rated(pitch,integral_error,dt)

   implicit none

   real, intent(inout) :: pitch
   real, intent(inout) :: integral_error

   real, intent(in) :: dt

   real :: dpitch
   real :: dpitch_max


! Reset integral state when leaving the above-rated region.
   integral_error = 0.0


! Return pitch gradually toward pitch_min.
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
!      Q = sum[(r x F) . e_axis]       [LB torque]
!
!      Qphys = Q * torque_conv         [N m]
!
!      Paero = Qphys * omega           [W]
!
! Paero is stored and used by the pitch controller at the next
! controller update.
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
       size(power_initialized) /= nturb) then

      write(*,*) 'ERROR turbine_controller_feedback: controller arrays'
      write(*,*) 'nturb                   = ',nturb
      write(*,*) 'size(paero_last)         = ',size(paero_last)
      write(*,*) 'size(power_initialized)  = ',size(power_initialized)
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
!      Qphys = QLB * rho_c * L_c^5 / T_c^2
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
! Physical torque [N m].
!-----------------------------------------------------------------------
      torque_phys = torque*torque_conv


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
          allocated(pitch_integral) .and. &
          allocated(power_initialized)) then

         if (size(paero_last)       == nturb .and. &
             size(pitch_integral)   == nturb .and. &
             size(power_initialized)== nturb) return

      endif

   endif


!-----------------------------------------------------------------------
! Allocate/reallocate controller state.
!-----------------------------------------------------------------------
   if (allocated(paero_last))        deallocate(paero_last)
   if (allocated(pitch_integral))    deallocate(pitch_integral)
   if (allocated(power_initialized)) deallocate(power_initialized)

   allocate(paero_last(nturb))
   allocate(pitch_integral(nturb))
   allocate(power_initialized(nturb))

   paero_last        = 0.0
   pitch_integral    = 0.0
   power_initialized = .false.

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
