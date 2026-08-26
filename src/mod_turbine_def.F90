module mod_turbine_def
   implicit none

!------------------------------------------------------------
! Blade geometry definitions
!
! hubradius   : Radial position of the inner edge of the blade,
!               measured from the rotor centre [m].
!
! rotorradius : Aerodynamic blade length from the blade root to
!               the blade tip [m].
! radius      : hubradius+rotorradius
!
!               Total rotor radius = hubradius + rotorradius
! ratedrpm      = 12.1       ! [rev/min]
! ratedpower    = 5.0e6      ! [W]
! pitch_min     = 0.0        ! [deg]
! pitch_max     = 30.0       ! [deg]
! pitchrate_max = 8.0        ! [deg/s]
!
! DC          : Spanwise length of the blade element [m].
!
! RELM        : Radial location of the centre of the blade element,
!               measured from the rotor centre [m].
!
!               RELM(i) = hubradius
!                         + SUM(DC(1:i-1))
!                         + DC(i)/2
!
! Chord       : Local airfoil chord length [m], i.e. the distance
!               from the leading edge to the trailing edge.
!
! Twist       : Local aerodynamic twist angle [deg], defined about
!               the blade spanwise (+x) axis. Positive twist follows
!               the conventional aerodynamic definition and decreases
!               the local angle of attack according to
!
!                   alpha = phi - twist - pitch
!
!               where phi is the local flow angle and pitch is the
!               blade pitch angle.
!
!               If geometry is imported from SIMA, whose aerodynamic
!               twist convention has the opposite sign, the imported
!               twist values must first be negated.
!
! foilfile(i) : File containing the aerodynamic coefficient table
!               associated with blade element i.
!
! Each blade element has its own foil table in foil_h(:,i). This may
! duplicate tables when several blade elements use the same airfoil,
! but makes the indexing between blade elements and foil data direct.
!------------------------------------------------------------

   integer, parameter :: nblades = 3
   real,    save :: rotorradius    ! [m]
   real,    save :: hubradius      ! [m]
   real,    save :: radius         ! [m]
   real,    save :: radiusnd       ! [LU]
   integer, save :: iradius        ! [-]
   integer, save :: nrchords       ! number of blade elements
   integer, save :: nrd            ! maximum number of entries in a foil table
   real,    save :: ratedrpm       ! [rev/min]
   real,    save :: ratedpower     ! [W]
   real,    save :: pitch_min      ! [deg]
   real,    save :: pitch_max      ! [deg]
   real,    save :: pitchrate_max  ! [deg/s]

   real, allocatable, save :: dc(:)
   real, allocatable, save :: relm(:)
   real, allocatable, save :: chord(:)
   real, allocatable, save :: twist(:)

   character(len=200), allocatable, save :: foilfile(:)

   ! ------------------------------
   ! Foil tables
   ! ------------------------------

   type foildata
      real :: deg
      real :: cl
      real :: cd
      real :: cm
   end type foildata

   ! Host copies (read from files)
   integer,        allocatable, save :: nrang_h(:)
   type(foildata), allocatable, save :: foil_h(:,:)

   ! Device copies (used by GPU kernels)
   integer,        allocatable, save :: nrang(:)
   type(foildata), allocatable, save :: foil(:,:)

#ifdef _CUDA
   attributes(device) :: nrang, foil
#endif

contains

!=======================================================================
! Read turbine blade geometry and allocate turbine arrays.
!
! Expected input format (NREL5MW):
!
! rotorradius   = 61.500
! hubradius     = 1.500
! radius        = 63.0
! nrchords      = 17
! nrd           = 150
! ratedrpm      = 12.1       ! [rev/min]
! ratedpower    = 5.0e6      ! [W]
! pitch_min     = 0.0        ! [deg]
! pitch_max     = 30.0       ! [deg]
! pitchrate_max = 8.0        ! [deg/s]
!
! Foil          DC      RELM   Chord      Twist   foilfile
! foil_1    2.73330   2.86670  3.54200  13.3080   Cylinder1.dat
! ...
!=======================================================================
subroutine turbine_def(turbname)

   implicit none

   character(len=*), intent(in) :: turbname

   integer :: iu, ios, i
   character(len=256) :: line
   character(len=32)  :: foilname
   character(len=200)  :: filename

   real :: rinner, router
   real, parameter :: tol = 1.0e-3

   filename=trim(turbname)//'/turbine_def.in'
   !------------------------------------------------------------
   ! Open turbine blade definition
   !------------------------------------------------------------
   open(newunit=iu, file=trim(filename), status='old', action='read', iostat=ios)

   if (ios /= 0) then
      write(*,*) 'ERROR: unable to open turbine file: ', trim(filename)
      stop
   end if

   !------------------------------------------------------------
   ! Read scalar parameters
   !------------------------------------------------------------

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_real_parameter(line, rotorradius)

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_real_parameter(line, hubradius)

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_integer_parameter(line, nrchords)

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_integer_parameter(line, nrd)

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_real_parameter(line, ratedrpm)

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_real_parameter(line, ratedpower)

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_real_parameter(line, pitch_min)

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_real_parameter(line, pitch_max)

   read(iu,'(A)',iostat=ios) line
   if (ios /= 0) call turbine_read_error(filename)
   call read_real_parameter(line, pitchrate_max)

   !------------------------------------------------------------
   ! Basic checks
   !------------------------------------------------------------
   if (nrchords <= 0) then
      write(*,*) 'ERROR: nrchords must be positive'
      stop
   end if

   if (nrd <= 0) then
      write(*,*) 'ERROR: nrdata must be positive'
      stop
   end if

   if (hubradius < 0.0 .or. rotorradius <= 0.0) then
      write(*,*) 'ERROR: invalid turbine radius specification'
      stop
   end if

   !------------------------------------------------------------
   ! Allocate blade geometry arrays
   !------------------------------------------------------------
   allocate(dc(nrchords))
   allocate(relm(nrchords))
   allocate(chord(nrchords))
   allocate(twist(nrchords))
   allocate(foilfile(nrchords))

   !------------------------------------------------------------
   ! Allocate airfoil arrays.
   !
   ! There is one table for each blade element.
   !------------------------------------------------------------
   allocate(nrang_h(nrchords))
   allocate(nrang(nrchords))

   allocate(foil_h(nrd,nrchords))
   allocate(foil(nrd,nrchords))

   nrang_h = 0

   !------------------------------------------------------------
   ! Find table header
   !------------------------------------------------------------
   do
      read(iu,'(A)',iostat=ios) line

      if (ios /= 0) then
         write(*,*) 'ERROR: turbine blade table not found'
         stop
      end if

      if (index(line,'Foil') > 0 .and. &
          index(line,'Chord') > 0) exit
   end do

   !------------------------------------------------------------
   ! Read blade elements
   !
   ! Example:
   !
   ! foil_1  2.73330  2.86670  3.54200  13.3080  Cylinder1.dat
   !------------------------------------------------------------
   do i = 1, nrchords

      read(iu,*,iostat=ios) foilname, dc(i), relm(i), &
                            chord(i), twist(i), foilfile(i)

      if (ios /= 0) then
         write(*,*) 'ERROR reading blade element ', i
         stop
      end if

      if (dc(i) <= 0.0) then
         write(*,*) 'ERROR: DC must be positive for blade element ', i
         stop
      end if

      if (chord(i) <= 0.0) then
         write(*,*) 'ERROR: chord must be positive for blade element ', i
         stop
      end if

   end do

   close(iu)

   !------------------------------------------------------------
   ! Geometry consistency checks
   !------------------------------------------------------------

   ! Inner edge of first element
   rinner = relm(1) - 0.5*dc(1)

   ! Outer edge of last element
   router = relm(nrchords) + 0.5*dc(nrchords)

   ! Expected blade tip radius
   radius = hubradius + rotorradius

   if (abs(rinner - hubradius) > tol) then
      write(*,*) 'WARNING: first blade element does not start at hubradius'
      write(*,*) '         hubradius        = ', hubradius
      write(*,*) '         first inner edge = ', rinner
   end if

   if (abs(router - radius) > tol) then
      write(*,*) 'WARNING: blade elements do not end at rotor tip'
      write(*,*) '         expected tip = ', radius
      write(*,*) '         table tip    = ', router
   end if

   !------------------------------------------------------------
   ! Print summary
   !------------------------------------------------------------
   write(*,*)
   write(*,*) 'Turbine blade definition read from: ', trim(filename)
   write(*,*) 'Blade length       = ', rotorradius,             '[m]'
   write(*,*) 'Hub radius         = ', hubradius,               '[m]'
   write(*,*) 'Total rotor radius = ', hubradius + rotorradius, '[m]'
   write(*,*) 'Number of elements = ', nrchords
   write(*,*) 'Maximum foil data  = ', nrd
   write(*,*) 'Rated RPM          = ', ratedrpm,                '[rev/min]'
   write(*,*) 'Rated power        = ', ratedpower,              '[W]'
   write(*,*) 'Minimum pitch      = ', pitch_min
   write(*,*) 'Maximum pitch      = ', pitch_max
   write(*,*) 'Maximum pitchrate  = ', pitchrate_max

   write(*,*)
   write(*,'(A)') &
      '   i       dc       relm      chord      twist    foil file'

   do i = 1, nrchords
      write(*,'(I4,4F11.4,2X,A)') &
           i, dc(i), relm(i), chord(i), twist(i), trim(foilfile(i))
   end do

end subroutine turbine_def


!=======================================================================
! Error reading turbine definition file
!=======================================================================
subroutine turbine_read_error(filename)

   implicit none

   character(len=*), intent(in) :: filename

   write(*,*) 'ERROR reading turbine definition file: ', trim(filename)
   stop

end subroutine turbine_read_error


!=======================================================================
! Extract a real value appearing after "="
!
! Example:
!
!    rotorradius = 61.500
!=======================================================================
subroutine read_real_parameter(line, value)

   implicit none

   character(len=*), intent(in) :: line
   real,             intent(out) :: value

   integer :: ipos, ios

   ipos = index(line,'=')

   if (ipos == 0) then
      write(*,*) 'ERROR: expected "=" in line: ', trim(line)
      stop
   end if

   read(line(ipos+1:),*,iostat=ios) value

   if (ios /= 0) then
      write(*,*) 'ERROR reading parameter from: ', trim(line)
      stop
   end if

end subroutine read_real_parameter


!=======================================================================
! Extract an integer value appearing after "="
!=======================================================================
subroutine read_integer_parameter(line, value)

   implicit none

   character(len=*), intent(in) :: line
   integer,          intent(out) :: value

   integer :: ipos, ios

   ipos = index(line,'=')

   if (ipos == 0) then
      write(*,*) 'ERROR: expected "=" in line: ', trim(line)
      stop
   end if

   read(line(ipos+1:),*,iostat=ios) value

   if (ios /= 0) then
      write(*,*) 'ERROR reading parameter from: ', trim(line)
      stop
   end if

end subroutine read_integer_parameter

subroutine turbine_read_foils()

   implicit none

   integer :: ic

   do ic = 1, nrchords
      write(*,*) 'Reading foil ', ic, ': ', trim(foilfile(ic))
      call turbine_read_foil_file(trim(foilfile(ic)), ic)
   end do

   nrang = nrang_h
   foil  = foil_h

end subroutine turbine_read_foils

subroutine turbine_read_foil_file(foilfile, ic)

   use m_readinfile, only : turbname
   implicit none

   character(len=*), intent(in) :: foilfile
   integer,          intent(in) :: ic

   character(len=100) filename

   integer :: iu, ios, iosdata, n
   real :: alpha, cl, cd, cm
   character(len=512) :: line

   filename=trim(turbname)//'/'//trim(foilfile)
   open(newunit=iu, file=trim(filename), status='old', action='read', iostat=ios)

   if (ios /= 0) then
      write(*,*) 'ERROR opening foil file: ', trim(filename)
      stop
   end if

   n = 0

   do

      read(iu,'(A)',iostat=ios) line

      if (ios < 0) exit

      if (ios > 0) then
         write(*,*) 'ERROR reading foil file: ', trim(filename)
         stop
      end if

      line = adjustl(line)

      ! Skip blank lines
      if (len_trim(line) == 0) cycle

      ! Skip comment/header lines
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      ! Read aerodynamic data
      read(line,*,iostat=iosdata) alpha, cl, cd, cm

      if (iosdata /= 0) then
         write(*,*) 'ERROR reading aerodynamic data in: ', trim(filename)
         write(*,*) 'Line: ', trim(line)
         stop
      end if

      n = n + 1

      if (n > nrd) then
         write(*,*) 'ERROR: too many entries in foil file: ', trim(filename)
         write(*,*) 'Maximum allowed = ', nrd
         stop
      end if

      foil_h(n,ic)%deg = alpha
      foil_h(n,ic)%cl  = cl
      foil_h(n,ic)%cd  = cd
      foil_h(n,ic)%cm  = cm

   end do

   close(iu)

   if (n < 2) then
      write(*,*) 'ERROR: too few aerodynamic data points in: ', &
                 trim(filename)
      stop
   end if

   nrang_h(ic) = n
   write(*,'(A,A,A,I4)') 'Read ', trim(filename), ': ', n

end subroutine turbine_read_foil_file
end module mod_turbine_def
