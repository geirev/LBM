module m_readinfile
! Simulation parameters
   integer  nt0      ! First timestep
   integer  nt1      ! Last timestep
   integer  iout     ! number of steps between outputs 0, 0+iout, ...
   integer  ifout    ! number of steps between outputs 0, 0+iout, ...
   integer  ibnd     ! Type of bondary condition in i direction
   integer  jbnd     ! Type of bondary condition in i direction
   integer  kbnd     ! Type of bondary condition in k direction
   logical  lpseudo  ! Add smooth pseudorandom peturbations to initial rho
   real     uini     ! Initial u-velocity
   real     rho0     ! Average density
   real     rhoa     ! Imposed density gradient for ibnd=2 case
   real     tau      ! Collision timescale 0.6
!  character(len=20) :: experiment='airfoil'
   character(len=20) :: experiment='cylinder'
!  character(len=20) :: experiment='sphere'
!  character(len=20) :: experiment='cube'
!  character(len=20) :: experiment='disks'

   type physconv
      real rho
      real length
      real time
      real vel
      real dynvisc
      real kinvisc
      real eddvisc
   end type
   real machnr
   type(physconv) p2l

contains
subroutine readinfile
implicit none

   character(len=2) ca
   character(len=3) :: version='1.0'
   character(len=3) ver
   logical ex
   real reynoldsnr
   real nu
   real newtau
   real gridrn

! reading input data
   inquire(file='infile.in',exist=ex)
   if (.not.ex) then
      print '(a)','Did not find inputfile infile.in...'
      stop
   endif
   print '(a)','--------------------------------------------------------------------------------'
   open(10,file='infile.in')
      read(10,'(tr1,a)')ver
      if (version /= ver) then
         print *,'Update version of infile.in to:',version,ver
         stop
      endif
      read(10,*)experiment         ; print '(a,a)',       'experiment        = ',trim(experiment)
      read(10,*)nt0                ; print '(a,i8)',      'nt0               = ',nt0
      read(10,*)nt1                ; print '(a,i8)',      'nt1               = ',nt1
      read(10,*)iout               ; print '(a,i8)',      'iout              = ',iout
      read(10,*)ifout              ; print '(a,i8)',      'ifout             = ',ifout
      read(10,*)ibnd               ; print '(a,i8)',      'ibnd              = ',ibnd
      read(10,*)jbnd               ; print '(a,i8)',      'jbnd              = ',jbnd
      read(10,*)kbnd               ; print '(a,i8)',      'kbnd              = ',kbnd
      read(10,*)uini               ; print '(a,f8.3)',    'uini (ini u vel)  = ',uini
      read(10,*)rho0               ; print '(a,f8.3)',    'rho0 (latt dens)  = ',rho0
      read(10,*)rhoa               ; print '(a,f8.3)',    'rhoa (pres grad)  = ',rhoa
      read(10,'(1x,l1)')lpseudo    ; print '(a,tr7,l1)',  'lpseudo           = ',lpseudo
      read(10,*)tau                ; print '(a,f8.3)',    'tau               = ',tau
      read(10,*)p2l%rho            ; print '(a,f8.3,a)',  'air density       = ',p2l%rho,    ' [kg/m^3]'   ! 1.225 is Air density at 15C and  101.325 kPa  (kg/m^3)
      read(10,*)p2l%length         ; print '(a,f8.3,a)',  'grid cell size    = ',p2l%length, ' [m]'
      read(10,*)p2l%vel            ; print '(a,f8.3,a)',  'wind velocity     = ',p2l%vel,    ' [m/s]'
      read(10,*)p2l%eddvisc        ; print '(a,f8.3,a)',  'eddy viscosity    = ',p2l%eddvisc,' [m^2/s]'
   close(10)

   print *
!  Compute time step from lattice size and velocity
   p2l%time=p2l%length/p2l%vel
   print '(a,f8.3,a)',  'p2l_time         = ',p2l%time     ,' [s]'

!  Compupte Reynolds number from rotor of 20 lattice cells
   reynoldsnr=(20.0*p2l%length)*p2l%vel/p2l%eddvisc
   print '(a,i8,a)',    'Reynolds num     = ',nint(reynoldsnr)   ,' [ ]'


!  Comupte Eddy viscosity from input tau
   nu=(1.0/3.0)*(tau - 0.5) * (p2l%length**2/p2l%time)           !(7.14)
   print '(a,f8.3,a)',  'nu from tau      = ',nu           ,' [m^2/s]'

   newtau=0.5+ p2l%eddvisc * (p2l%time/p2l%length**2) / 3.0
   print '(a,f8.3,a)', 'tau from eddy visc= ',newtau       ,' [ ]'
   print '(a,f8.3)',   'constraint:    tau> ',0.5+0.3/8.0

   nu=(1.0/3.0)*(newtau - 0.5) * (p2l%length**2/p2l%time)
   print '(a,f8.3,a)', 'nu new            = ',nu           ,' [m^2/s]'

   gridrn= p2l%length*p2l%vel/p2l%eddvisc
   print '(a,f8.3,a)', 'Grid Reynolds num = ',gridrn       ,' [ ]'

   print *
   print '(a,g12.4)','Error terms:'
   print '(a,g12.4)','Spatial discretization errors proportional to dx^2       :', p2l%length**2
   print '(a,g12.4)','Time    discretization errors proportional to dt^2       :', p2l%time**2
   print '(a,g12.4)','Compressibility        errors proportional to dt^2/dx**2 :', p2l%time**2/p2l%length**2
   print '(a,g12.4)','BGK truncation         errors proportional to (tau-0.5)*2:', (newtau-0.5)**2

end subroutine
end module

