module mod_dimensions

!windfarm big D=32
!   integer, parameter :: nx = 928          ! 29D resolution x-dir (east)
!   integer, parameter :: ny = 121          ! 3D  resolution y-dir (north)
!   integer, parameter :: nz = 121          ! 3D  resolution z-dir (up)

!windfarm
!  integer, parameter :: nx = 928          ! 928 resolution x-dir (east)
!  integer, parameter :: ny = 96          ! resolution y-dir (north)
!  integer, parameter :: nz = 96          ! resolution z-dir (up)

!city
!   integer, parameter :: nx = 200          ! 928 resolution x-dir (east)
!   integer, parameter :: ny = 96           ! resolution y-dir (north)
!   integer, parameter :: nz = 96            ! resolution z-dir (up)

!minicity
    integer, parameter :: nx = 200          ! 928 resolution x-dir (east)
    integer, parameter :: ny = 120          ! resolution y-dir (north)
    integer, parameter :: nz = 1            ! resolution z-dir (up)

!cylinder and airfoil
!  integer, parameter :: nx = 400          ! resolution x-dir (east)
!  integer, parameter :: ny = 100          ! resolution y-dir (north)
!  integer, parameter :: nz = 5            ! resolution z-dir (up)

! mini            
! integer, parameter :: nx = 3          ! resolution x-dir (east)
! integer, parameter :: ny = 3          ! resolution y-dir (north)
! integer, parameter :: nz = 3          !  resolution z-dir (up)


end module
