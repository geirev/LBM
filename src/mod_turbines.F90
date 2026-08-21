module mod_turbines
! Global turbine types, parameters and shared arrays

   use mod_dimensions, only : nx, ny, nz, nyg

#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : mpi_rank, mpi_nprocs, j_start, j_end
#endif

   implicit none

   real, parameter :: pi  = acos(-1.0)
   real, parameter :: pi2 = 2.0*pi


!------------------------------------------------------------
! Stores turbine-specific configuration information.
!
! Blade geometry (relm, dc, chord, twist) and foil tables are
! stored globally in mod_turbine_def and are shared by all
! turbines.
!------------------------------------------------------------
   type turbine_t

      integer :: imodel           ! 0 = actuator line (ALM)
                                  ! 1 = rotating actuator disk (ADM-R)

      real    :: xhub             ! hub x-position (global grid index)
      real    :: yhub             ! hub y-position (global grid index)
      real    :: zhub             ! hub z-position (global grid index)

      real    :: radius           ! total rotor radius, non-dimensional
                                  ! = (hubradius + rotorradius)/p2l%length

      integer :: iradius          ! rotor radius in number of grid cells

      integer :: nblades          ! number of blades

      real    :: theta            ! rotor azimuth (rad)
      real    :: yaw              ! yaw angle (rad)
      real    :: tilt             ! tilt angle (rad)

      real    :: omegand          ! non-dim angular speed, omega*p2l%time
      real    :: pitchangle       ! collective pitch (deg)
      integer :: tiploss          ! tip-loss flag

      integer :: nazim            ! number of azimuthal samples per radius,
                                  ! ADM-R only

   end type turbine_t


!------------------------------------------------------------
! One actuator sample point per chord per blade per turbine
! in global coordinates.
!------------------------------------------------------------
   type :: point_t

      integer :: iturb            ! turbine index
      integer :: iblade           ! blade/sample index
      integer :: ichord           ! blade element index

      real    :: xg, yg, zg       ! global coordinates (grid units)

      real    :: dc               ! local radial segment length
      real    :: chord            ! local chord width
      real    :: relm             ! local radial position

      integer :: foil             ! foil-table index; normally = ichord

      real    :: twist            ! local aerodynamic twist (deg)

      real    :: yaw, tilt        ! yaw/tilt (rad)
      real    :: theta            ! rotor azimuth (rad)
      real    :: pitch            ! collective pitch (deg)
      real    :: omegand          ! non-dim angular speed

      real    :: force_scale      ! multiplies Fvec after
                                  ! turbine_compute_bladeforce;
                                  ! 1.0 for ALM,
                                  ! nblades/nazim for ADM-R

   end type point_t


!------------------------------------------------------------
! Global storage
!------------------------------------------------------------
   type(turbine_t), allocatable :: turbines(:)
   type(point_t),   allocatable :: points_global(:)


! Point block limits (not used yet)
   integer :: t_imin, t_imax
   integer :: t_jmin, t_jmax
   integer :: t_kmin, t_kmax


! Time-filtered local upstream velocities for each turbine
   real,    allocatable, save :: uavg_f(:), vavg_f(:), wavg_f(:)
   logical, allocatable, save :: windfilter_initialized(:)

end module mod_turbines
