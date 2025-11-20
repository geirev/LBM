module mod_turbines
! Global turbine types, parameters and shared arrays
   use mod_dimensions, only : nx, ny, nz, nyg
#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : mpi_rank, mpi_nprocs, j_start, j_end
#endif
   use mod_nrel5mw,   only : nrchords, relm, dc, chord, nfoil, twist, hubradius, rotorradius
   implicit none
   real, parameter :: pi  = 3.14159265358979323846
   real, parameter :: pi2 = 2.0*pi

!   Stores NREL-5MW style configuration information.
   type turbine_t
      real    :: xhub             ! hub x-position (global grid index)
      real    :: yhub             ! hub y-position (global grid index)
      real    :: zhub             ! hub z-position (global grid index)

      real    :: radius           ! rotor radius (non-dimensional, /p2l%length)
      integer :: iradius          ! rotor radius in number of grid cells

      integer :: nblades          ! number of blades
      integer :: nchords          ! chords per blade

      real    :: relm(nrchords)   ! radial positions (non-dimensional)
      real    :: dc   (nrchords)  ! radial segment length (non-dimensional)
      real    :: chord(nrchords)  ! chord width (non-dimensional)
      real    :: twist(nrchords)  ! twist angle per chord (deg)
      integer :: nfoil(nrchords)  ! foil index

      real    :: theta            ! rotor azimuth (rad)
      real    :: yaw              ! yaw angle (rad)
      real    :: tilt             ! tilt angle (rad)

      real    :: omegand          ! non-dim angular speed, omega*p2l%time
      real    :: pitchangle       ! collective pitch (deg)
      integer :: tiploss          ! tip-loss flag
   end type turbine_t

!   One actuator sample point per chord per blade per turbine in global coordinates.
   type :: point_t
      integer :: iturb            ! turbine index
      integer :: iblade           ! blade index
      integer :: ichord           ! chord index

      real    :: xg, yg, zg       ! global coordinates (in grid units)
      real    :: dc               ! local radial segment
      real    :: chord            ! local chord width
      real    :: relm             ! local radial position
      integer :: foil             ! foil index
      real    :: twist            ! local twist (deg)
      real    :: yaw, tilt        ! yaw/tilt (rad)
      real    :: theta            ! rotor azimuth (rad)
      real    :: pitch            ! pitch (deg)
      real    :: omegand          ! non-dim angular speed
   end type point_t

! Global storage (used by high-level drivers if desired)
   type(turbine_t), allocatable :: turbines(:)
   type(point_t),   allocatable :: points_global(:)

! Storage for forcing on lattice nodes
   real                         :: F_turb(3,0:nx+1,0:ny+1,0:nz+1)

! point block limits (not used yet)
   integer t_imin,t_imax,t_jmin,t_jmax,t_kmin,t_kmax
end module mod_turbines

