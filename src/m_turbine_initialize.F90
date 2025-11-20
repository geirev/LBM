!==============================================================
!  m_turbine_initalize.F90
!  Turbine initialization: read config and fill turbines(:)
!==============================================================
module m_turbine_initialize
   use mod_turbines
   use m_readinfile, only : nturbines, p2l, ipos, jpos, kpos, yaw, tilt, turbrpm, pitchangle, tipspeedratio, itiploss
   use m_nrelreadfoil
   implicit none
contains

!--------------------------------------------------------------
!  subroutine turbine_initialize
!
!  PURPOSE:
!    Allocate and initialize the global array turbines(:)
!    using the configuration from m_readinfile and NREL-5MW
!    geometry from mod_nrel5mw.
!
!  CALL:
!    call turbine_initialize()
!
!  SIDE EFFECTS:
!    - Allocates mod_turbines::turbines(:)
!    - Calls nrelreadfoil() to load airfoil tables
!--------------------------------------------------------------
subroutine turbine_initialize()
   implicit none
   integer :: n
   real    :: radius, omega

   if (.not. allocated(turbines))  allocate(turbines(nturbines))

   ! Rotation speed in radians / second
   omega = pi2 * turbrpm / 60.0

   do n = 1, nturbines
      ! Hub position in global grid indices
      turbines(n)%xhub = real(ipos(n))
      turbines(n)%yhub = real(jpos(n))
      turbines(n)%zhub = real(kpos(n))

      ! Rotor geometry (non-dimensionalized by p2l%length)
      radius              = rotorradius + hubradius
      turbines(n)%radius  = radius / p2l%length
      turbines(n)%iradius = nint(turbines(n)%radius)

!      print '(a,f8.2,i4)', 'Rotor radius=   ', turbines(n)%radius,  turbines(n)%iradius
!      print '(a,f8.2)',    'Rotor diameter=', 2.0 * turbines(n)%radius

      ! Blade discretization
      turbines(n)%nblades = 3
      turbines(n)%nchords = nrchords

      turbines(n)%relm(1:)  = relm(1:)   / p2l%length
      turbines(n)%dc  (1:)  = dc  (1:)   / p2l%length
      turbines(n)%chord(1:) = chord(1:)  / p2l%length
      turbines(n)%twist(1:) = twist(1:)
      turbines(n)%nfoil(1:) = nfoil(1:)

      ! Orientation & dynamics
      turbines(n)%theta      = 0.0
      turbines(n)%yaw        = (yaw(n)/360.0)*pi2
      turbines(n)%tilt       = (tilt(n)/360.0)*pi2
      turbines(n)%omegand    = omega * p2l%time
      turbines(n)%pitchangle = pitchangle
      turbines(n)%tiploss    = itiploss
   end do


   ! Load airfoil tables for all foils used
   call nrelreadfoil()
end subroutine turbine_initialize

end module m_turbine_initialize
