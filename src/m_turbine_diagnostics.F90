!=======================================================================
! Turbine aerodynamic diagnostics
!
! Computes for each turbine:
!
!   - total axial thrust
!   - total rotor torque
!   - torque distribution over blade elements
!   - mechanical rotor power
!   - Ct, Cq and Cp based on local upstream velocity
!   - Ct, Cq and Cp based on prescribed reference inflow velocity
!   - actual tip-speed ratio
!   - angular velocity and physical RPM
!
!-----------------------------------------------------------------------
! LATTICE AND PHYSICAL UNITS
!-----------------------------------------------------------------------
!
! The turbine calculation is performed internally in lattice units.
!
! Physical conversion factors are:
!
!      p2l%rho       density scale           [kg/m^3]
!      p2l%length    lattice-cell length     [m]
!      p2l%time      lattice timestep        [s]
!      p2l%vel       velocity scale          [m/s]
!
! with
!
!      p2l%time = p2l%length/p2l%vel
!
! The lattice variables therefore have the following dimensions:
!
!      position, radius, relm     [LU]
!      velocity                   [LU/TS]
!      density                    [rho_LB]
!      omegand                    [rad/TS]
!      force                      [rho_LB LU^4/TS^2]
!      torque                     [rho_LB LU^5/TS^2]
!      power                      [rho_LB LU^5/TS^3]
!
! The physical conversions are:
!
!      U [m/s] =
!           U_LB * p2l%vel
!
!      R [m] =
!           R_LB * p2l%length
!
!      omega [rad/s] =
!           omegand / p2l%time
!
!      F [N] =
!           F_LB * p2l%rho*p2l%length^4/p2l%time^2
!
!      Q [N m] =
!           Q_LB * p2l%rho*p2l%length^5/p2l%time^2
!
!      P [W] =
!           P_LB * p2l%rho*p2l%length^5/p2l%time^3
!
!-----------------------------------------------------------------------
! FORCE AND TORQUE CONVENTIONS
!-----------------------------------------------------------------------
!
! Fvec_global contains the aerodynamic force exerted by the flow on
! the turbine blade.
!
! The rotor basis is right-handed:
!
!      e_axis x e_rot = e_tan
!
! Consequently
!
!      Q = sum[(r x F) . e_axis]
!
! and independently
!
!      Q = sum[r F_tan]
!
! where positive torque acts in the direction of positive rotor
! rotation.
!
! Mechanical turbine power is therefore
!
!      P = Q*omega
!
! and positive P corresponds to power extracted by the rotor.
!
! The dimensionless coefficients are
!
!      Ct = T / (0.5*rho*U^2*A)
!
!      Cq = Q / (0.5*rho*U^2*A*R)
!
!      Cp = P / (0.5*rho*U^3*A)
!
! with
!
!      lambda = omega*R/U
!
! giving the exact consistency relation
!
!      Cp = lambda*Cq
!
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

   real, parameter :: pi      = acos(-1.0)
   real, parameter :: pi2     = 2.0*pi
   real, parameter :: rho_ref = 1.0

   real :: e_axis(3), e1(3), e2(3)
   real :: e_tan(3)

   real :: rvec(3)
   real :: fpoint(3)
   real :: torque_vec(3)

   real :: thrust
   real :: torque
   real :: torque_rft
   real :: tq

   real :: power

   real :: ct
   real :: cq
   real :: cp

   real :: ct_ref
   real :: cq_ref
   real :: cp_ref

   real :: lambda_actual

   real :: area
   real :: radius
   real :: radius_phys

   real :: ulocal
   real :: ulocal_phys

   real :: omega
   real :: rpm

   real :: force_conv
   real :: torque_conv
   real :: power_conv

   real :: thrust_phys
   real :: torque_phys
   real :: power_phys

   real, allocatable :: torque_chord(:)

!-----------------------------------------------------------------------
! Physical conversion factors.
!-----------------------------------------------------------------------
   force_conv  = p2l%rho*p2l%length**4/p2l%time**2
   torque_conv = p2l%rho*p2l%length**5/p2l%time**2
   power_conv  = p2l%rho*p2l%length**5/p2l%time**3


   allocate(torque_chord(nrchords))


   do n = 1, nturbines

      thrust       = 0.0
      torque       = 0.0
      torque_rft   = 0.0
      torque_chord = 0.0


!-----------------------------------------------------------------------
! Rotor basis.
!-----------------------------------------------------------------------
      call turbine_rotor_basis(turbines_in(n)%yaw, &
                               turbines_in(n)%tilt, &
                               e_axis, e1, e2)


!-----------------------------------------------------------------------
! Integrate actuator-point forces.
!-----------------------------------------------------------------------
      do p = 1, np

         if (points_global(p)%iturb /= n) cycle

         fpoint(:) = Fvec_global(:,p)

         rvec(1) = points_global(p)%xg - turbines_in(n)%xhub
         rvec(2) = points_global(p)%yg - turbines_in(n)%yhub
         rvec(3) = points_global(p)%zg - turbines_in(n)%zhub

         e_tan(:) = -sin(points_global(p)%theta)*e1(:) + &
                      cos(points_global(p)%theta)*e2(:)

         ! Axial rotor thrust [LB force].
         thrust = thrust + dot_product(fpoint,e_axis)

         ! Rotor torque from Q = (r x F) . e_axis [LB torque].
         torque_vec(1) = rvec(2)*fpoint(3) - rvec(3)*fpoint(2)
         torque_vec(2) = rvec(3)*fpoint(1) - rvec(1)*fpoint(3)
         torque_vec(3) = rvec(1)*fpoint(2) - rvec(2)*fpoint(1)

         tq = dot_product(torque_vec,e_axis)

         torque = torque + tq

         ! Independent torque calculation Q = r F_tan [LB torque].
         torque_rft = torque_rft + points_global(p)%relm * &
                      dot_product(fpoint,e_tan)

         ic = points_global(p)%ichord
         torque_chord(ic) = torque_chord(ic) + tq

      enddo


!-----------------------------------------------------------------------
! Torque consistency checks.
!-----------------------------------------------------------------------
      if (abs(torque-torque_rft) > &
          1.0e-5*max(1.0,abs(torque))) then

         write(*,'(A,ES14.6,2X,A,ES14.6,2X,A,ES12.4)') &
              'Torque(r x F) [LB]=', torque, &
              'Torque(r Ft) [LB]=',  torque_rft, &
              'difference [LB]=',    torque-torque_rft

         write(*,*) 'ERROR: torque diagnostic mismatch'
         stop

      endif

      if (abs(torque-sum(torque_chord)) > &
          1.0e-5*max(1.0,abs(torque))) then

         write(*,*) 'ERROR: radial torque diagnostic mismatch'
         write(*,*) 'torque [LB]           = ', torque
         write(*,*) 'sum torque_chord [LB] = ', sum(torque_chord)

         stop

      endif


!-----------------------------------------------------------------------
! Performance coefficients.
!-----------------------------------------------------------------------
      if (windfilter_initialized(n)) then

         ulocal = sqrt(uavg_f(n)**2 + &
                       vavg_f(n)**2 + &
                       wavg_f(n)**2)

         ulocal_phys = ulocal*p2l%vel

         radius      = turbines_in(n)%radius
         radius_phys = radius*p2l%length

         area = pi*radius*radius


         if (ulocal > 1.0e-12) then

            ! Actual tip-speed ratio [-].
            lambda_actual = turbines_in(n)%omegand*radius/ulocal

            ! Mechanical rotor power [LB power].
            power = torque*turbines_in(n)%omegand

            !---------------------------------------------------------
            ! Coefficients based on filtered local upstream wind.
            !---------------------------------------------------------
            ct = thrust / (0.5*rho_ref*ulocal**2*area)

            cq = torque / (0.5*rho_ref*ulocal**2*area*radius)

            cp = power  / (0.5*rho_ref*ulocal**3*area)

            !---------------------------------------------------------
            ! Coefficients based on prescribed reference inflow.
            !---------------------------------------------------------
            ct_ref = thrust / (0.5*rho_ref*uini**2*area)

            cq_ref = torque / (0.5*rho_ref*uini**2*area*radius)

            cp_ref = power  / (0.5*rho_ref*uini**3*area)


            !---------------------------------------------------------
            ! Consistency check:
            !
            !       Cp = lambda*Cq
            !---------------------------------------------------------
            if (abs(cp-lambda_actual*cq) > &
                1.0e-5*max(1.0,abs(cp))) then

               write(*,*) 'ERROR: Cp /= lambda*Cq'
               write(*,*) 'Cp        = ',cp
               write(*,*) 'lambda*Cq = ',lambda_actual*cq
               stop

            endif

         else

            lambda_actual = 0.0
            power         = 0.0

            ct     = 0.0
            cq     = 0.0
            cp     = 0.0

            ct_ref = 0.0
            cq_ref = 0.0
            cp_ref = 0.0

         endif


!-----------------------------------------------------------------------
! Convert rotor dynamics and aerodynamic loads to physical units.
!-----------------------------------------------------------------------
         omega = turbines_in(n)%omegand/p2l%time
         rpm   = omega*60.0/pi2

         thrust_phys = thrust*force_conv
         torque_phys = torque*torque_conv
         power_phys  = power*power_conv


!-----------------------------------------------------------------------
! Diagnostic output.
!
! U       : filtered local upstream wind speed        [m/s]
! R       : physical rotor radius                     [m]
! omega   : angular velocity                          [rad/s]
! RPM     : rotor speed                               [rev/min]
! TSR     : actual tip-speed ratio                    [-]
! TSR*    : prescribed target tip-speed ratio         [-]
! Ct      : local-wind thrust coefficient             [-]
! Cq      : local-wind torque coefficient             [-]
! Cp      : local-wind power coefficient              [-]
! Ctr     : reference-inflow thrust coefficient       [-]
! Cqr     : reference-inflow torque coefficient       [-]
! Cpr     : reference-inflow power coefficient        [-]
! T       : axial rotor thrust                        [kN]
! Q       : rotor torque                              [kN m]
! P       : mechanical rotor power                    [MW]
!-----------------------------------------------------------------------

         write(*,'(A,I3, &
                   &2X,A,F7.2,A, &
                   &2X,A,F7.2,A, &
                   &2X,A,F8.3,A, &
                   &2X,A,F7.2,A, &
                   &2X,A,F6.2,A, &
                   &2X,A,F6.2,A, &
                   &2X,A,F7.3,A, &
                   &2X,A,F7.3,A, &
                   &2X,A,F7.3,A, &
                   &2X,A,F7.3,A, &
                   &2X,A,F7.3,A, &
                   &2X,A,F7.3,A, &
                   &2X,A,F9.2,A, &
                   &2X,A,F10.2,A, &
                   &2X,A,F9.3,A)') &
              'Turbine ', n, &
              'U=',       ulocal_phys,          ' m/s', &
              'R=',       radius_phys,          ' m', &
              'omega=',   omega,                ' rad/s', &
              'RPM=',     rpm,                  ' rev/min', &
              'TSR=',     lambda_actual,        ' [-]', &
              'TSR*=',    tipspeedratio,        ' [-]', &
              'Ct=',      ct,                   ' [-]', &
              'Cq=',      cq,                   ' [-]', &
              'Cp=',      cp,                   ' [-]', &
              'Ctr=',     ct_ref,               ' [-]', &
              'Cqr=',     cq_ref,               ' [-]', &
              'Cpr=',     cp_ref,               ' [-]', &
              'T=',       thrust_phys/1.0e3,    ' kN', &
              'Q=',       torque_phys/1.0e3,    ' kN m', &
              'P=',       power_phys/1.0e6,     ' MW'

      endif


!-----------------------------------------------------------------------
! Radial torque distribution.
!-----------------------------------------------------------------------
      cycle

      write(*,'(A)') &
           '   ic      r/R [-]        Q [kN m]       fraction [-]'

      do ic = 1, nrchords

         write(*,'(I5,F12.4,ES16.6,F16.6)') &
              ic, &
              relm(ic)/turbines_in(n)%radius, &
              torque_chord(ic)*torque_conv/1.0e3, &
              torque_chord(ic)/torque

      enddo

   enddo


   deallocate(torque_chord)

end subroutine turbine_diagnostics

end module m_turbine_diagnostics
