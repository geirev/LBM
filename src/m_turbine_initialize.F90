!==============================================================
!  m_turbine_initalize.F90
!  Turbine initialization: read config and fill turbines(:)
!==============================================================
module m_turbine_initialize
   use mod_turbines
   use mod_turbine_def
   use m_readinfile, only : nturbines, p2l, xpos,ypos,zpos,yaw,tilt,turbinename,turbrpm, pitchangle, itiploss, alm_adm, nazim
   implicit none
contains

!--------------------------------------------------------------
!  subroutine turbine_initialize
!
!  PURPOSE:
!    Allocate and initialize the global array turbines(:)
!    using the configuration from m_readinfile
!
!  CALL:
!    call turbine_initialize()
!
!  SIDE EFFECTS:
!    - Allocates mod_turbines::turbines(:)
!    - Calls read_foil() to load airfoil tables
!--------------------------------------------------------------
subroutine turbine_initialize()
   implicit none
   integer :: n
   real    :: radius, omega

! Upstream tubine velocity filtering
   allocate(uavg_f(nturbines))
   allocate(vavg_f(nturbines))
   allocate(wavg_f(nturbines))
   allocate(windfilter_initialized(nturbines))
   uavg_f = 0.0
   vavg_f = 0.0
   wavg_f = 0.0
   windfilter_initialized = .false.


   if (.not. allocated(turbines))  allocate(turbines(nturbines))

   ! Rotation speed in radians / second
   omega = pi2 * turbrpm / 60.0

   ! Blade geometry is read in meters; convert to lattice/grid units
   relm(1:)  = relm(1:)   / p2l%length
   dc  (1:)  = dc  (1:)   / p2l%length
   chord(1:) = chord(1:)  / p2l%length

   do n = 1, nturbines

      turbines(n)%name = trim(turbinename(n))
      turbines(n)%imodel = alm_adm

      ! Hub position in global coordinates (meter)
      turbines(n)%xhub = xpos(n)/p2l%length
      turbines(n)%yhub = ypos(n)/p2l%length
      turbines(n)%zhub = zpos(n)/p2l%length

      ! Rotor geometry (non-dimensionalized by p2l%length)
      radius              = rotorradius + hubradius
      turbines(n)%radius  = radius / p2l%length
      turbines(n)%iradius = nint(turbines(n)%radius)

!      print '(a,f8.2,i4)', 'Rotor radius=   ', turbines(n)%radius,  turbines(n)%iradius
!      print '(a,f8.2)',    'Rotor diameter=', 2.0 * turbines(n)%radius

      ! Blade discretization
      turbines(n)%nblades = 3
      turbines(n)%nrchords = nrchords


      ! Orientation & dynamics
      turbines(n)%theta      = 0.0
      turbines(n)%yaw        = (yaw(n)/360.0)*pi2
      turbines(n)%tilt       = (tilt(n)/360.0)*pi2
      turbines(n)%omegand    = omega * p2l%time
      turbines(n)%pitchangle = pitchangle
      turbines(n)%tiploss    = itiploss

      ! Actuator disk samples per radius
      turbines(n)%nazim = nazim
   end do

   ! Load aerodynamic coefficient tables
   call turbine_read_foils()

end subroutine turbine_initialize
end module m_turbine_initialize
