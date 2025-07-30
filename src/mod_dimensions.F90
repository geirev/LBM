module mod_dimensions
!windfarm
!   integer, parameter :: nx = 928          ! 928 resolution x-dir (east)
!   integer, parameter :: ny = 96           ! resolution y-dir (north)
!   integer, parameter :: nz = 96           ! resolution z-dir (up)

!city
   integer, parameter :: nx = 200          ! 928 resolution x-dir (east)
   integer, parameter :: ny = 96           ! resolution y-dir (north)
   integer, parameter :: nz = 96            ! resolution z-dir (up)

!cylinder and airfoil
!  integer, parameter :: nx = 400          ! 928 resolution x-dir (east)
!  integer, parameter :: ny = 100          ! resolution y-dir (north)
!  integer, parameter :: nz = 5            ! resolution z-dir (up)

! Numbers of threads per block are read from infile.in
   integer :: ntx                          ! Number of threads per block in x-direction
   integer :: nty                          ! Number of threads per block in y-direction
   integer :: ntz                          ! Number of threads per block in z-direction

! 42.9 43.16 seconds
!  integer, parameter :: ntx = 16
!  integer, parameter :: nty = 4
!  integer, parameter :: ntz = 8

! 43.0
!  integer, parameter :: ntx = 16
!  integer, parameter :: nty = 2
!  integer, parameter :: ntz = 8

! 43.6 seconds
!  integer, parameter :: ntx = 16
!  integer, parameter :: nty = 4
!  integer, parameter :: ntz = 6

! 45.4 seconds
!   integer, parameter :: ntx = 8
!   integer, parameter :: nty = 4
!   integer, parameter :: ntz = 8

! 45.9 seconds
!  integer, parameter :: ntx = 12
!  integer, parameter :: nty = 4
!  integer, parameter :: ntz = 8

! 45.68 seconds
!  integer, parameter :: ntx = 4
!  integer, parameter :: nty = 8
!  integer, parameter :: ntz = 16

! 59 seconds
!  integer, parameter :: ntx = 8
!  integer, parameter :: nty = 16
!  integer, parameter :: ntz = 4

!  52.7 seconds
!  integer, parameter :: ntx = 32
!  integer, parameter :: nty = 4
!  integer, parameter :: ntz = 4


!  53.1 seconds
!   integer, parameter :: ntx = 8
!   integer, parameter :: nty = 8
!   integer, parameter :: ntz = 4

! 56.3 seconds
!   integer, parameter :: ntx = 16
!   integer, parameter :: nty = 8
!   integer, parameter :: ntz = 4

!
!  48.8 seconds
!   integer, parameter :: ntx = 32
!   integer, parameter :: nty = 4
!   integer, parameter :: ntz = 1

!  50.4 seconds
!   integer, parameter :: ntx = 16
!   integer, parameter :: nty = 4
!   integer, parameter :: ntz = 4
!
!  45.2 seconds
!   integer, parameter :: ntx = 8
!   integer, parameter :: nty = 8
!   integer, parameter :: ntz = 8

!  49.2 seconds
!    integer, parameter :: ntx = 16
!    integer, parameter :: nty = 4
!    integer, parameter :: ntz = 2



   integer, parameter :: sz=8                ! size of a real
   integer, parameter :: nl = 27             ! number of componets of D3Q27
   integer, parameter :: ieps = 5            ! number of gridcells for smoothing actuatorline forcing in i-dir
   integer, parameter :: nrturb = 1000       ! number of precomputed batches if inflow turbulence for u,v,w, and rho
end module
