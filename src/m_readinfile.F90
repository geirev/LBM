module m_readinfile

   use mod_dimensions, only : nx, nyg, nz

   implicit none

!===============================================================================
!  Simulation parameters
!===============================================================================

   integer  nt0             ! First timestep
   integer  nt1             ! Last timestep
   integer  iout            ! Number of steps between outputs: 0, iout, ...
   integer  iprt1           ! Start of high-frequency output interval
   integer  iprt2           ! End of high-frequency output interval
   integer  dprt            ! High-frequency output interval
   logical  ltesting        ! Write minimal Tecplot files if true
   integer  irestart        ! Number of steps between restart files
   integer  itecout         ! Tecplot output format:
                            !   0 = full files
                            !   2 = solution variables only

!===============================================================================
!  Boundary conditions
!===============================================================================

   integer  ibnd            ! Boundary-condition type in i-direction
   integer  jbnd            ! Boundary-condition type in j-direction
   integer  kbnd            ! Boundary-condition type in k-direction

   logical  lsponge         ! Add sponge layer along open boundaries
   real     tau_max         ! Maximum tau in sponge layer
   integer  nsponge_i       ! Sponge-layer thickness in i-direction
   integer  nsponge_j       ! Sponge-layer thickness in j-direction

!===============================================================================
!  Inflow
!===============================================================================

   logical  inflowturbulence   ! Add smooth pseudo-random inflow perturbations
   integer  nrturb             ! Number of precomputed turbulence batches
   real     timecor            ! Time correlation of turbulence forcing fields
   real     turb_length        ! Length scales of turbulence forcing fields
   real     turbulence_ampl    ! Amplitude of inflow turbulence

   real     uini               ! Initial absolute velocity
   real     udir               ! Initial velocity direction [degrees]
   real     zref               ! Reference/hub height [m]
   real     z0                 ! Surface roughness length [m]
   logical  wind_profile       ! Apply vertical wind profile
   logical  wall_model         ! Apply a wall model on the lower boundary

!===============================================================================
!  Physical and numerical parameters
!===============================================================================

   real     rho0               ! Reference density
   real     tauin              ! Collision timescale
   real     kinevisc           ! Kinematic viscosity
   real     machnr             ! Mach number

   character(len=20) experiment   ! Experiment name

   type physconv
      real rho
      real length
      real time
      real vel
      real visc
   end type physconv

   type(physconv) p2l


!===============================================================================
!  GPU configuration
!===============================================================================
   integer  ntx              ! Threads per block in x-direction
   integer  nty              ! Threads per block in y-direction
   integer  ntz              ! Threads per block in z-direction

!===============================================================================
!  Averaging and diagnostics
!===============================================================================

   logical  laveraging       ! Compute full-grid averages (memory demanding)
   integer  avestart         ! Iteration at which averaging starts
   integer  avesave          ! Iteration at which averages are saved

   logical  ldump            ! Dump diagnostic files to disk
   logical  lmeasurements    ! Used in data-assimilation experiments
   logical  ltiming          ! Enable kernel timing; false avoids synchronizations

!===============================================================================
!  Turbines
!===============================================================================

   character(len=20) turbname       ! Turbine definition identifier
   integer           nturbines      ! Number of turbines
   integer           actuator_model ! 0 = actuator line, 1 = actuator disk
   integer           nazim          ! Number of azimuthal points per radius
   integer           localwind      ! 1 = local wind, 2 = external udir,
                                    ! otherwise use infile yaw

   real     filter_time       ! Time filter for local upstream wind [s]
   real     dtcontrol         ! Yaw-controller update interval [s]
   real     yawrate_max       ! Maximum yaw rate [deg/s]
   real     yaw_deadband      ! Yaw-controller deadband [deg]
   real     pitchangle        ! Imposed pitch angle [deg]
   real     turbrpm           ! Imposed turbine speed [rpm]
   real     tipspeedratio     ! Imposed tip-speed ratio
   real     powerloss         ! Turbine drivetrain/generator power-loss factor

   real, allocatable :: xpos(:), ypos(:), zpos(:)
   real, allocatable :: yaw(:), tilt(:)
   character(len=10), allocatable :: turbinename(:)

!===============================================================================
!  Collision, SGS and atmospheric boundary layer
!===============================================================================

   integer  ihrr             ! 1 = regularized R(fneq) scheme
   integer  ibgk             ! 2/3 = second/third-order BGK equilibrium
   integer  ivreman          ! 1 = Vreman SGS model
   integer  iablvisc         ! 0 = none, 1 = mechanical, 2 = buoyancy
   real     ablheight        ! Atmospheric boundary-layer height [m]
   integer  istable          ! -1 = unstable, 0 = neutral, 1 = stable
   real     smagorinsky      ! Smagorinsky constant
   integer  iturb


contains


subroutine readinfile()

   use m_mkinfile

   implicit none

   character(len=1) :: ver
   logical          :: ex
   real             :: gridrn
   integer          :: n


!===============================================================================
!  Check execution directory
!===============================================================================

   inquire(file='main.F90',exist=ex)

   if (ex) then
      stop 'You are executing boltzmann in the src/ directory'
   endif


!===============================================================================
!  Check for input file
!===============================================================================

   inquire(file='infile.in',exist=ex)

   if (.not.ex) then

      print '(a)','Did not find input file infile.in'
      print '(a)','Generating new template infile.in.....'

      call mkinfile()

      print '(a)','Please edit infile.in and relaunch boltzmann'

      stop

   endif


!===============================================================================
!  Read input file
!===============================================================================

   print '(a)', &
      '--------------------------------------------------------------------------------'

   open(10,file='infile.in')

!-------------------------------------------------------------------------------
!  Diagnostics and GPU configuration
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,*,err=100) ltiming
   print '(a,tr7,l1)', &
      'ltiming           = ',ltiming

   read(10,*,err=100) ltesting
   print '(a,tr7,l1)', &
      'ltesting          = ',ltesting

   read(10,*,err=100) ldump
   print '(a,tr7,l1)', &
      'ldump             = ',ldump

   read(10,*,err=100) ntx,nty,ntz
   print '(a,3i4)', &
      'threads per block = ',ntx,nty,ntz

   read(10,*,err=100) lmeasurements
   print '(a,tr7,l1)', &
      'lmeasurements     = ',lmeasurements


!-------------------------------------------------------------------------------
!  Numerical model
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,*,err=100) experiment
   print '(a,a)', &
      'experiment        = ',trim(experiment)

   read(10,*,err=100) ibgk
   print '(a,i1)', &
      'BGK order of feq  = ',ibgk

   read(10,*,err=100) ihrr
   print '(a,i1)', &
      'HRR regularization= ',ihrr

   read(10,*,err=100) ivreman,smagorinsky
   print '(a,i1,a,f10.4)', &
      'Vreman mixing     = ',ivreman, &
      ' Smagorinsky=',smagorinsky


!-------------------------------------------------------------------------------
!  Atmospheric boundary layer
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,*,err=100) iablvisc
   print '(a,i1)', &
      'ABL: iablvisc     = ',iablvisc

   read(10,*,err=100) ablheight
   print '(a,f10.4)', &
      'ABL: ablheight    = ',ablheight

   read(10,*,err=100) istable
   print '(a,i3)', &
      'ABL: stability    = ',istable


!-------------------------------------------------------------------------------
!  Timestepping and output
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,*,err=100) nt0
   print '(a,i8)', &
      'nt0               = ',nt0

   read(10,*,err=100) nt1
   print '(a,i8)', &
      'nt1               = ',nt1

   if (nt1 <= nt0) then
      stop 'readinfile: nt1 <= nt0'
   endif

   read(10,*,err=100) iout
   print '(a,i8)', &
      'iout              = ',iout

   read(10,*,err=100) irestart
   print '(a,i8)', &
      'irestart          = ',irestart

   read(10,*,err=100) iprt1,iprt2,dprt
   print '(a,3i8)', &
      'iprt1, iprt2, dprt= ',iprt1,iprt2,dprt

   read(10,*,err=100) itecout
   print '(a,i8)', &
      'itecout           = ',itecout


!-------------------------------------------------------------------------------
!  Boundary conditions
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,*,err=100) ibnd
   print '(a,i8)', &
      'ibnd              = ',ibnd

   read(10,*,err=100) jbnd
   print '(a,i8)', &
      'jbnd              = ',jbnd

   read(10,*,err=100) kbnd
   print '(a,i8)', &
      'kbnd              = ',kbnd

   lsponge   = .false.
   tau_max   = 0.7
   nsponge_i = 15
   nsponge_j = 15


!-------------------------------------------------------------------------------
!  Inflow
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,*,err=100) uini,udir
   print '(a,2(f8.3,a))', &
      'inflow (uini,udir)= ',uini,' [m/s]',udir,' [degrees]'

   read(10,*,err=100) wind_profile,wall_model,zref,z0
   print '(a,tr7,2l1,tr2,f10.2,f10.6)', &
      'wind profile, wallmodel, zref,z0 = ',wind_profile,wall_model,zref,z0

   timecor=0.98
   turb_length=30.0       ! [m]
   read(10,*,err=100) inflowturbulence,turbulence_ampl,nrturb
   print '(a,tr7,l1,tr2,g13.5,tr2,i5)', &
      'inflowturbulence  = ', &
      inflowturbulence,turbulence_ampl,nrturb


!-------------------------------------------------------------------------------
!  Physical scaling
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,*,err=100) kinevisc
   print '(a,f8.3,a)', &
      'Kinematic viscos  = ',kinevisc,' [m^2/s]'

   read(10,*,err=100) p2l%rho
   print '(a,f8.3,a)', &
      'air density       = ',p2l%rho,' [kg/m^3]'

   read(10,*,err=100) p2l%length
   print '(a,f8.3,a)', &
      'grid cell size    = ',p2l%length,' [m]'

   read(10,*,err=100) p2l%vel
   print '(a,f8.3,a)', &
      'wind velocity     = ',p2l%vel,' [m/s]'

   ! Convert dimensional inflow velocity to lattice units.
   uini = uini/p2l%vel

   print '(a,f8.3,a)', &
      'Non-dim uinflow   = ',uini,' [] Should be less than 0.2'


!-------------------------------------------------------------------------------
!  Averaging
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,'(1x,l1)',err=100) laveraging
   print '(a,tr7,l1)', &
      'laveraging        = ',laveraging

   read(10,*,err=100) avestart
   print '(a,i8)', &
      'avestart iteration= ',avestart

   read(10,*,err=100) avesave
   print '(a,i8)', &
      'avesave iteration = ',avesave


!-------------------------------------------------------------------------------
!  Turbines
!-------------------------------------------------------------------------------

   read(10,'(a)',err=100) ver

   read(10,*,err=100) nturbines
   print '(a,i8)', &
      'Num of turbines   = ',nturbines

   if (nturbines > 0) then

      allocate( &
         xpos(nturbines), &
         ypos(nturbines), &
         zpos(nturbines), &
         yaw(nturbines),  &
         tilt(nturbines), &
         turbinename(nturbines))

      read(10,*,err=100) turbname
      print '(a,a)', &
         'Turbine type      = ',trim(turbname)

      read(10,*,err=100) actuator_model
      print '(a,i8)', &
         'Actuator model    = ',actuator_model

      read(10,*,err=100) nazim
      print '(a,i8)', &
         'nazim in ADM      = ',nazim

      read(10,*,err=100) localwind
      print '(a,i8)', &
         'localwind         = ',localwind

      read(10,*,err=100) dtcontrol
      print '(a,f8.3,a)', &
         'dtcontrol         = ',dtcontrol,' [s]'

      read(10,*,err=100) yawrate_max
      print '(a,f8.3,a)', &
         'yawrate_max       = ',yawrate_max,' [deg/s]'

      read(10,*,err=100) yaw_deadband
      print '(a,f8.3,a)', &
         'yaw_deadband      = ',yaw_deadband,' [deg]'

      read(10,*,err=100) filter_time
      print '(a,f8.3,a)', &
         'filter_time       = ',filter_time,' [s]'

      read(10,*,err=100) pitchangle
      print '(a,f8.3,a)', &
         'Pitch angle       = ',pitchangle,' [deg]'

      read(10,*,err=100) turbrpm
      print '(a,f8.3,a)', &
         'RPM for act.line  = ',turbrpm,' [rotations/min]'

      read(10,*,err=100) tipspeedratio
      print '(a,f8.3,a)', &
         'Tipspeed ratio    = ',tipspeedratio,' []'

      read(10,*,err=100) powerloss
      print '(a,f8.3,a)', &
         'Powerloss factor  = ',powerloss,' []'

      read(10,'(a)',err=100) ver

      do n=1,nturbines

         read(10,*,err=100) &
            iturb,turbinename(n), &
            xpos(n),ypos(n),zpos(n),yaw(n),tilt(n)

         print '(a,x,i4,x,3a,x,5f10.2)', &
            'Turbine',iturb, &
            '(',trim(turbinename(n)),'):', &
            xpos(n),ypos(n),zpos(n),yaw(n),tilt(n)

         if (xpos(n) < 0.0 .or. &
             xpos(n) > real(nx )*p2l%length .or. &
             ypos(n) < 0.0 .or. &
             ypos(n) > real(nyg)*p2l%length .or. &
             zpos(n) < 0.0 .or. &
             zpos(n) > real(nz )*p2l%length) then

            write(*,*)
            write(*,*) 'ERROR: turbine hub lies outside model domain'

            write(*,'(a,i5,2x,a)') &
               'Turbine : ',n,trim(turbinename(n))

            write(*,'(a,3f12.3,a)') &
               'Hub     : ',xpos(n),ypos(n),zpos(n),' [m]'

            write(*,'(a,2f12.3,a)') &
               'x-range : ',0.0,real(nx)*p2l%length,' [m]'

            write(*,'(a,2f12.3,a)') &
               'y-range : ',0.0,real(nyg)*p2l%length,' [m]'

            write(*,'(a,2f12.3,a)') &
               'z-range : ',0.0,real(nz)*p2l%length,' [m]'

            error stop

         endif

      enddo

   else

      print '(a)','Running without wind turbines'

   endif

   close(10)


!===============================================================================
!  Conversion factors
!===============================================================================

   print *
   print '(a)','Conversion factors:'

   print '(a,f12.6,a)', &
      'p2l%length = ',p2l%length,' [m]'

   p2l%time = p2l%length/p2l%vel

   print '(a,f12.6,a)', &
      'p2l%time   = ',p2l%time,' [s]'

   print '(a,f12.6,a)', &
      'p2l%vel    = ',p2l%vel,' [m/s]'

   print '(a,f12.6,a)', &
      'p2l%rho    = ',p2l%rho,' [kg/m^3]'

   p2l%visc = p2l%length**2/p2l%time

   print '(a,f12.4,a)', &
      'p2l%visc   = ',p2l%visc,' [m^2/s]'

   print *


!===============================================================================
!  Compute nondimensional tau from dimensional kinematic viscosity
!===============================================================================

   print '(a,g13.6,a)', &
      'Kinematic visc       = ',kinevisc,' [m^2/s]'

   tauin = 0.5 + 3.0*kinevisc/p2l%visc

   print '(a,g13.6,a)', &
      'tau from kinevisc    = ',tauin,' [ ]'


!===============================================================================
!  Nondimensional kinematic viscosity used in fequil
!===============================================================================

   kinevisc = kinevisc/p2l%visc

   print '(a,g13.6,a)', &
      'Non-dim kinevisc     = ',kinevisc,' [ ]'


!===============================================================================
!  Grid-cell Reynolds number
!===============================================================================

   gridrn = p2l%length*uini*p2l%vel / &
            (kinevisc*p2l%visc)

   print '(a,i12,a)', &
      'cell-Reynolds num    = ',nint(gridrn),' [ ]'


!===============================================================================
!  Mach number
!===============================================================================

   print '(a,f8.3,a)', &
      'Mach number (u/c)    = ',uini*p2l%vel/330.0,' [ ]'


!===============================================================================
!  Error estimates
!===============================================================================

   print *
   print '(a)','Error terms:'

   print '(a,g12.4)', &
      'Spatial discretization errors proportional to dx^2       :', &
      p2l%length**2

   print '(a,g12.4)', &
      'Time    discretization errors proportional to dt^2       :', &
      p2l%time**2

   print '(a,g12.4)', &
      'Compressibility        errors proportional to dt^2/dx**2 :', &
      p2l%time**2/p2l%length**2

   print '(a,g12.4)', &
      'BGK truncation       errors proportional to (tauin-0.5)^2:', &
      (tauin-0.5)**2


!===============================================================================
!  Reference density
!===============================================================================

   rho0 = 1.0

   return


!===============================================================================
!  Input-file error
!===============================================================================

100 continue

   close(10)

   call system('mv -i infile.in infile_backup.in')
   call mkinfile()

   print *
   print '(a)','ERROR while reading infile.in'
   print '(a)','The original input file has been moved to infile_backup.in'
   print '(a)','A new template infile.in has been generated.'
   print '(a)','Please check the input file and relaunch boltzmann.'

   stop


end subroutine readinfile

end module m_readinfile
