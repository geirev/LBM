module m_readinfile
! Simulation parameters
   integer  nt0            ! First timestep
   integer  nt1            ! Last timestep
   integer  iout           ! number of steps between outputs 0, 0+iout, ...
   integer  iprt           ! Output every time steps of it <= iprt
   integer  ifout          ! number of steps between outputs 0, 0+iout, ...
   integer  ibnd           ! Type of bondary condition in i direction
   integer  jbnd           ! Type of bondary condition in i direction
   integer  kbnd           ! Type of bondary condition in k direction
   logical  lpseudo        ! Add smooth pseudorandom peturbations to initial rho
   logical  runexp         ! Add smooth pseudorandom peturbations to initial rho
   character(len=3) coll   ! Collision operator
   integer  ihrr           ! Option for collisions using HRR scheme
   real     uini           ! Initial u-velocity
   real     rho0           ! Average density
   real     rhoa           ! Imposed density gradient for ibnd=2 case
   real     tau            ! Collision timescale 0.6
   character(len=20) experiment ! experiment name

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

   integer nturbines  ! Number of turbines in model
   real turbrad       ! Turbine radius in meters
   integer radii      ! Turbine radius in grid cells
   real turbblock     ! The wind velocity reduction induced by a turbine
   integer, allocatable ::  ipos(:),jpos(:),kpos(:) ! Turbine locations

contains
subroutine readinfile
   implicit none

   character(len=3) ca
   character(len=3) :: version='1.0'
   character(len=3) ver
   logical ex
   real reynoldsnr
   real nu
   real newtau
   real gridrn
   integer n

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
      read(10,'(1x,l1)')runexp     ; print '(a,tr7,l1)',  'runexp            = ',runexp
      read(10,*)experiment         ; print '(a,a)',       'experiment        = ',trim(experiment)
      read(10,*)coll,ihrr          ; print '(3a,i1)',     'Collision operator= ',coll,' with ihrr option: ',ihrr
      if ((coll /= 'HRR') .and. (coll /= 'BGK')) stop 'invalid collision operator'
      read(10,*)nt0                ; print '(a,i8)',      'nt0               = ',nt0
      read(10,*)nt1                ; print '(a,i8)',      'nt1               = ',nt1
      read(10,*)iout               ; print '(a,i8)',      'iout              = ',iout
      read(10,*)iprt               ; print '(a,i8)',      'iprt              = ',iprt
      read(10,*)ifout              ; print '(a,i8)',      'ifout             = ',ifout
      read(10,*)ibnd               ; print '(a,i8)',      'ibnd              = ',ibnd
      read(10,*)jbnd               ; print '(a,i8)',      'jbnd              = ',jbnd
      read(10,*)kbnd               ; print '(a,i8)',      'kbnd              = ',kbnd
      read(10,*)uini               ; print '(a,f8.3,a)',  'uini (inflow uvel)= ',uini,       ' [m/s]'
      read(10,*)rho0               ; print '(a,f8.3,a)',  'rho0 (latt dens)  = ',rho0,       ' []'
      read(10,*)rhoa               ; print '(a,f8.3,a)',  'rhoa (pres grad)  = ',rhoa,       ' []'
      read(10,'(1x,l1)')lpseudo    ; print '(a,tr7,l1)',  'lpseudo           = ',lpseudo
      read(10,*)tau                ; print '(2(a,f8.3))', 'tau               = ',tau,        ' []  constraint: > ',0.5+0.3/8.0
      read(10,*)p2l%rho            ; print '(a,f8.3,a)',  'air density       = ',p2l%rho,    ' [kg/m^3]'   ! 1.225 is Air density at 15C and  101.325 kPa  (kg/m^3)
      read(10,*)p2l%length         ; print '(a,f8.3,a)',  'grid cell size    = ',p2l%length, ' [m]'
      read(10,*)p2l%vel            ; print '(a,f8.3,a)',  'wind velocity     = ',p2l%vel,    ' [m/s]'
      uini=uini/p2l%vel            ; print '(a,f8.3,a)',  'Non-dim uinflow   = ',uini,       ' [] Should be less that 0.2'

      read(10,'(a)')ca
      if (ca /= '#-T') then
         print *,'#-T: error in infile.in'
         stop
      endif
      print *
      read(10,*)nturbines              ; print '(a,i8)',          'Num of turbines   = ',nturbines
      if (nturbines > 0) then
         allocate(ipos(nturbines), jpos(nturbines), kpos(nturbines))
         read(10,*)turbrad             ; print '(a,f8.3,a)',      'turbine radius    = ',turbrad,  ' [m]'
         radii=nint(turbrad/p2l%length); print '(a,i8,a)',        'turbine radii     = ',radii,    ' [grid cells]'
         read(10,*)turbblock           ; print '(a,f8.3)',        'turb blockage     = ',turbblock
         do n=1,nturbines
            read(10,*)ipos(n)
            read(10,*)jpos(n)
            read(10,*)kpos(n)
            print '(a,i4,a,3i4)', '(ijk)-pos for turbine  = ',n,' : ',ipos(n),jpos(n),kpos(n)
         enddo
      else 
         print '(a)','Running without wind turbines'
      endif

   close(10)

   print *
   print '(a)','Conversion factors:'
   print '(a,f12.4,a)','C_l=  ',p2l%length,   ' [m]'
   p2l%time=p2l%length/p2l%vel
   print '(a,f12.4,a)','C_t=  ',p2l%time,     ' [s]'
   print '(a,f12.4,a)','C_u=  ',p2l%vel,      ' [m/s]'
   print '(a,f12.4,a)','C_rho=',p2l%rho,      ' [kg/m^3]'
   print *

!  Compute dimensional eddy viscosity from input non-dimensional tau
   p2l%eddvisc=(1.0/3.0)*(tau - 0.5) * (p2l%length**2/p2l%time)     !(7.14)
   print '(a,g13.6,a)',  'eddy visc from tau = ',p2l%eddvisc  ,' [m^2/s]'

!  Compupte Reynolds number from rotor of radii lattice cells
   reynoldsnr=(2.0*real(radii)*p2l%length)*uini*p2l%vel/p2l%eddvisc
   print '(a,i12,a)',    'Reynolds num       = ',nint(reynoldsnr)   ,' [ ]'

!  Compupte grid cell Reynolds number
   gridrn= p2l%length*uini*p2l%vel/p2l%eddvisc
   print '(a,i12,a)',    'cell-Reynolds num  = ',nint(gridrn)       ,' [ ]'

! Mach number
   print '(a,f8.3,a)',  'Mach number (u/c)  = ',uini*p2l%vel/330.0 ,' [ ]'

   print *
   print '(a,g12.4)','Error terms:'
   print '(a,g12.4)','Spatial discretization errors proportional to dx^2       :', p2l%length**2
   print '(a,g12.4)','Time    discretization errors proportional to dt^2       :', p2l%time**2
   print '(a,g12.4)','Compressibility        errors proportional to dt^2/dx**2 :', p2l%time**2/p2l%length**2
   print '(a,g12.4)','BGK truncation         errors proportional to (tau-0.5)*2:', (newtau-0.5)**2

   if (.not.runexp) stop 'runexp is false'

end subroutine
end module

