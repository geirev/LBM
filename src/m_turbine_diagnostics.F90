!=======================================================================
! Turbine aerodynamic diagnostics
!
! Computes for each turbine:
!
!   - total axial thrust
!   - total rotor torque
!   - torque distribution over blade elements
!   - power
!   - Ct and Cp based on local upstream velocity
!   - Ct and Cp based on prescribed reference inflow velocity
!   - actual tip-speed ratio and physical RPM
!
! Torque is evaluated independently in two ways:
!
!      Q = sum[(r x F) . e_axis]
!
! and
!
!      Q = -sum[r F_tan]
!
! providing a consistency check on the force and rotor-coordinate
! conventions.
!=======================================================================
module m_turbine_diagnostics
contains
subroutine turbine_diagnostics(turbines_in, points_global, Fvec_global, np)

   use mod_turbines, only : turbine_t, point_t, &
                            uavg_f, vavg_f, wavg_f, &
                            windfilter_initialized
   use mod_turbine_def, only : nrchords, relm
   use m_readinfile, only : nturbines, p2l, uini, tipspeedratio
   use m_turbine_rotor_basis

   implicit none

   type(turbine_t), intent(in) :: turbines_in(:)
   type(point_t),   intent(in) :: points_global(:)
   integer,         intent(in) :: np
   real,            intent(in) :: Fvec_global(3,np)

   integer :: n, p, ic

   real, parameter :: pi  = acos(-1.0)
   real, parameter :: pi2 = 2.0*pi
   real, parameter :: rho_ref = 1.0

   real :: e_axis(3), e1(3), e2(3)
   real :: e_tan(3)
   real :: rvec(3), fpoint(3), torque_vec(3)

   real :: thrust, torque, torque_rft, tq
   real :: power
   real :: ct, cp, ct_ref, cp_ref
   real :: lambda_actual
   real :: area, radius, ulocal
   real :: omega, rpm

   real, allocatable :: torque_chord(:)


   allocate(torque_chord(nrchords))

   do n = 1, nturbines

      thrust       = 0.0
      torque       = 0.0
      torque_rft   = 0.0
      torque_chord = 0.0

      ! Rotor basis is common to all actuator points of this turbine
      call turbine_rotor_basis(turbines_in(n)%yaw, &
                               turbines_in(n)%tilt, &
                               e_axis, e1, e2)

      !------------------------------------------------------------
      ! Integrate actuator-point forces
      !------------------------------------------------------------
      do p = 1, np

         if (points_global(p)%iturb /= n) cycle

         fpoint(:) = Fvec_global(:,p)

         ! Radius vector from hub to actuator point
         rvec(1) = points_global(p)%xg - turbines_in(n)%xhub
         rvec(2) = points_global(p)%yg - turbines_in(n)%yhub
         rvec(3) = points_global(p)%zg - turbines_in(n)%zhub

         ! Local tangential direction
         e_tan(:) = -sin(points_global(p)%theta)*e1(:) + &
                      cos(points_global(p)%theta)*e2(:)

         ! Axial thrust
         thrust = thrust + dot_product(fpoint,e_axis)

         ! Torque from (r x F) . e_axis
         torque_vec(1) = rvec(2)*fpoint(3) - rvec(3)*fpoint(2)
         torque_vec(2) = rvec(3)*fpoint(1) - rvec(1)*fpoint(3)
         torque_vec(3) = rvec(1)*fpoint(2) - rvec(2)*fpoint(1)

         tq = dot_product(torque_vec,e_axis)

         torque = torque + tq

         ! Independent torque check.
         ! Minus sign follows the present rotor-axis convention.
         torque_rft = torque_rft - points_global(p)%relm * &
                      dot_product(fpoint,e_tan)

         ! Torque contribution from this radial element
         ic = points_global(p)%ichord
         torque_chord(ic) = torque_chord(ic) + tq

      enddo


      !------------------------------------------------------------
      ! Torque consistency checks
      !------------------------------------------------------------
      if (abs(torque-sum(torque_chord)) > &
          1.0e-5*max(1.0,abs(torque))) then
         write(*,'(A,ES14.6,2X,A,ES14.6,2X,A,ES12.4)') &
           'Torque(r x F)=', torque, &
           'Torque(r Ft)=', torque_rft, &
           'difference=', torque-torque_rft

         write(*,*) 'ERROR: torque diagnostic mismatch'
         write(*,*) 'torque           = ', torque
         write(*,*) 'sum torque_chord = ', sum(torque_chord)
         stop

      endif


      !------------------------------------------------------------
      ! Performance coefficients
      !------------------------------------------------------------
      if (windfilter_initialized(n)) then

         ulocal = sqrt(uavg_f(n)**2 + &
                       vavg_f(n)**2 + &
                       wavg_f(n)**2)

         radius = turbines_in(n)%radius
         area   = pi*radius*radius

         if (ulocal > 1.0e-12) then

            lambda_actual = turbines_in(n)%omegand*radius/ulocal

            power = torque*turbines_in(n)%omegand

            ! Coefficients based on local upstream wind
            ct = thrust / (0.5*rho_ref*ulocal**2*area)
            cp = -power / (0.5*rho_ref*ulocal**3*area)

            ! Coefficients based on prescribed reference inflow
            ct_ref = thrust / (0.5*rho_ref*uini**2*area)
            cp_ref = -power / (0.5*rho_ref*uini**3*area)

         else

            lambda_actual = 0.0
            power  = 0.0
            ct     = 0.0
            cp     = 0.0
            ct_ref = 0.0
            cp_ref = 0.0

         endif

         omega = turbines_in(n)%omegand/p2l%time
         rpm   = omega*60.0/pi2

         write(*,'(A,I3,2X,A,F6.2,2X,A,F6.2,2X,A,F6.2,2X,A,F6.2, &
                  &2(2X,A,F6.3,2X,A,F6.3),2X,A,F10.3)') &
              'Turbine', n, &
              'U=',     ulocal*p2l%vel, &
              'RPM=',   rpm, &
              'TSR=',   lambda_actual, &
              'TSR*=',  tipspeedratio, &
              'Ct=',    ct, &
              'Cp=',    cp, &
              'Ctr=',   ct_ref, &
              'Cpr=',   cp_ref, &
              'power=', power

      endif


      !------------------------------------------------------------
      ! Radial torque distribution
      !------------------------------------------------------------
      cycle
      write(*,'(A)') '   ic      r/R        torque       fraction'

      do ic = 1, nrchords
         write(*,'(I5,F10.4,2F13.5)') &
              ic, &
              relm(ic)/turbines_in(n)%radius, &
              torque_chord(ic), &
              torque_chord(ic)/torque
      enddo

   enddo

   deallocate(torque_chord)

end subroutine turbine_diagnostics
end module m_turbine_diagnostics
