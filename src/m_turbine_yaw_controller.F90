module m_turbine_yaw_controller
   implicit none
contains

pure real function wrap_pi(angle)
   use mod_turbines, only : pi2
   implicit none
   real, intent(in) :: angle

   wrap_pi = modulo(angle + 0.5*pi2, pi2) - 0.5*pi2

end function wrap_pi


subroutine turbine_yaw_controller(dir, dtcontrol, nturb, yaw, yawrate_max, yaw_deadband)
   use mod_turbines, only : pi2
   implicit none

   integer, intent(in) :: nturb
   real, intent(in)    :: dir
   real, intent(in)    :: dtcontrol
   real, intent(in)    :: yawrate_max
   real, intent(in)    :: yaw_deadband
   real, intent(inout) :: yaw(nturb)

   integer :: n
   real    :: deg2rad
   real    :: target_yaw
   real    :: yaw_error
   real    :: dyaw
   real    :: yawrate_max_rad
   real    :: yaw_deadband_rad

   deg2rad = pi2/360.0

   ! dir, yawrate_max and yaw_deadband are supplied in degrees.
   ! Turbine yaw is stored in radians.
   target_yaw      = wrap_pi(dir*deg2rad)
   yawrate_max_rad = yawrate_max*deg2rad
   yaw_deadband_rad = yaw_deadband*deg2rad

   do n = 1, nturb

      yaw_error = wrap_pi(target_yaw - yaw(n))

      if (abs(yaw_error) > yaw_deadband_rad) then

         dyaw = sign(min(abs(yaw_error), &
                        yawrate_max_rad*dtcontrol), yaw_error)

         yaw(n) = wrap_pi(yaw(n) + dyaw)

      endif

   enddo

end subroutine turbine_yaw_controller

end module m_turbine_yaw_controller
