module mod_dimensions
! ny is the grid dimension in the northward direction.
! If the domain is tiled using MPI parallelization with domain decomposition in the y-direction
! ny is the tile dimension, while nyg is an ntiles*ny the global dimension in the y-direction.

!windfarm big D=32
!   integer, parameter :: nx = 110          !928           !110          !928          ! grid dimension x-dir (east)
!   integer, parameter :: ny = 120          ! Local-tile grid dimension y-dir (north)
!   integer, parameter :: ntiles = 4        ! Number of tiles in y direction
!   integer, parameter :: nyg = ntiles*ny   ! global grid dimension y-dir (north)
!   integer, parameter :: nz = 121          ! grid dimension z-dir (up)
!   integer, parameter :: ntracer = 0       ! Number of tracer fields (potential temperature etc)

!windturbine yaw test
   integer, parameter :: nx = 240          ! grid dimension x-dir (east)
   integer, parameter :: ny = 120          ! grid dimension y-dir (north)
   integer, parameter :: ntiles = 2        ! Number of tiles in y direction
   integer, parameter :: nyg = ntiles*ny   ! global grid dimension y-dir (north)
   integer, parameter :: nz = 60          ! grid dimension z-dir (up)
   integer, parameter :: ntracer = 0       ! Number of tracer fields (potential temperature etc)

!city
!  integer, parameter :: nx = 200          ! grid dimension x-dir (east)
!  integer, parameter :: ny = 96           ! grid dimension y-dir (north)
!  integer, parameter :: nyg = 96          ! grid dimension y-dir (north)
!  integer, parameter :: nz = 96           ! grid dimension z-dir (up)

!city2
!  integer, parameter :: nx = 200          ! grid dimension x-dir (east)
!  integer, parameter :: ny = 120          ! grid dimension y-dir (north)
!  integer, parameter :: ntiles = 1        ! Number of tiles in y direction
!  integer, parameter :: nyg = ntiles*ny   ! global grid dimension y-dir (north)
!  integer, parameter :: nz = 2            ! grid dimension z-dir (up)
!  integer, parameter :: ntracer = 0       ! Number of tracer fields (potential temperature etc)

!city3
! integer, parameter :: nx = 52           ! grid dimension x-dir (east)
! integer, parameter :: ny = 31           ! grid dimension y-dir (north)
! integer, parameter :: ntiles = 1        ! Number of tiles in y direction
! integer, parameter :: nyg = ntiles*ny   ! global grid dimension y-dir (north)
! integer, parameter :: nz = 20           ! grid dimension z-dir (up)
! integer, parameter :: ntracer = 0       ! Number of tracer fields (potential temperature etc)

!cylinder and airfoil
! integer, parameter :: nx = 400          !nx = 400  grid dimension x-dir (east)
! integer, parameter :: ny = 50           !ny = 100  grid dimension y-dir (north)
! integer, parameter :: nz = 5            !nz = 5    grid dimension z-dir (up)
! integer, parameter :: ntiles = 2        ! Number of tiles in y direction
! integer, parameter :: nyg = ntiles*ny   ! global grid dimension y-dir (north)
! integer, parameter :: ntracer = 0       ! Number of tracer fields (potential temperature etc)

! mini
! integer, parameter :: nx = 3            ! grid dimension x-dir (east)
! integer, parameter :: ny = 3            ! grid dimension y-dir (north)
! integer, parameter :: nyg =3            ! grid dimension y-dir (north)
! integer, parameter :: nz = 3            ! grid dimension z-dir (up)

!Barcelona_4m
! integer, parameter :: nx = 226           ! 225 grid dimension x-dir (east)
! integer, parameter :: ny = 218           ! 218 grid dimension y-dir (north)
! integer, parameter :: ntiles = 1         ! Number of tiles in y direction
! integer, parameter :: nyg = ntiles*ny    ! global grid dimension y-dir (north)
! integer, parameter :: nz = 33            ! 32 grid dimension z-dir (up)
! integer, parameter :: ntracer = 0        ! Number of tracer fields (potential temperature etc)

!Barcelona_2m
!  integer, parameter :: nx = 453           ! 225 grid dimension x-dir (east)
!  integer, parameter :: ny = 437           ! 218 grid dimension y-dir (north)
!  integer, parameter :: nz = 67            ! 32 grid dimension z-dir (up)
!  integer, parameter :: ntiles = 1         ! Number of tiles in y direction
!  integer, parameter :: nyg = ntiles*ny    ! global grid dimension y-dir (north)
!  integer, parameter :: ntracer = 0        ! Number of tracer fields (potential temperature etc)

end module
