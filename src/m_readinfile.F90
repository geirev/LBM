module m_readinfile
! Simulation parameters
   integer  nt0            ! First timestep
   integer  nt1            ! Last timestep
   integer  iout           ! number of steps between outputs 0, 0+iout, ...
   integer  iprt1          ! Output every dprt time steps of it <= iprt
   integer  iprt2          ! Output every dprt time steps of it <= iprt
   integer  dprt           ! delta high frequency output 
   logical  lprtmin        ! Print minimalistice plt file if true (no derived variables)
   integer  irestart       ! number of steps between restart files
   integer  ifout          ! number of steps between outputs 0, 0+iout, ...
   integer  ibnd           ! Type of bondary condition in i direction (ibnd=0 periodic, 1 inflow/outflow, 12
   integer  jbnd           ! Type of bondary condition in i direction
   integer  kbnd           ! Type of bondary condition in k direction
   logical  linipert       ! Add smooth pseudo-random peturbations to initial fields
   logical  inflowturbulence          ! Add smooth pseudo-random peturbations for turbulent inflow
   integer  nrturb         ! Number of precomputed batches of inflow turbulence for u,v,w, and rho
   real     turbulence_ampl! strength of inflow turbulence
   real     uini           ! Initial absolute velocity
   real     udir           ! Initial velocity direction in degrees
   real     rho0           ! Average density
   real     rhoa           ! Imposed density gradient for ibnd=2 case
   real     tauin          ! Collision timescale
   real     kinevisc       ! Kinematic viscosity (nondimensional used in fequil)
   character(len=20) experiment ! experiment name
   logical laveraging      ! Computes full averages ower the whole grid (memory demanding)
   logical laveturb        ! Computes  section averages with wind turbine according to Ashmut
   integer avestart        ! Iteration number for starting to compute averages
   integer avesave         ! Iteration number for ending averaging and saving averages
   integer itiploss        ! Tiploss(0-none, 1-Prandl, 2-Shen)
   integer :: ntx          ! Number of threads per block in x-direction
   integer :: nty          ! Number of threads per block in y-direction
   integer :: ntz          ! Number of threads per block in z-direction

   type physconv
      real rho
      real length
      real time
      real vel
      real visc
   end type
   real machnr
   type(physconv) p2l

   integer nturbines       ! Number of turbines in model
   real pitchangle         ! Imposed pitch angle
   real turbrpm            ! Imposed turbine RPM
   real tipspeedratio      ! Imposed tipspeed ratio
   integer, allocatable ::  ipos(:),jpos(:),kpos(:) ! Turbine locations
   integer  ihrr           ! Option (1) for regularized R(fneq) scheme
   integer  ibgk           ! Option (2,3) for second or third order BGK f^eq expansion
   integer  ivreman        ! Option (1) for subgridscale mixing using Vreman
   real smagorinsky        ! smagorinsky constant (0.15) used in subgridscale mixing
   integer iforce          ! Method for forcing scheme
                           !  iforce=1  !  Shan and Chen (1993)
                           !  iforce=8  !  Guo (2002)
                           !  iforce=10 !  Kupershtokh (2009)
                           !  iforce=12 !  Khazaeli et al. 2019

contains
subroutine readinfile()
   use mod_dimensions
   implicit none

   character(len=3) ca
   character(len=3) :: version='1.0'
   character(len=3) ver
   logical ex
   real gridrn
   integer n
   logical  runexp

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
      read(10,*)ntx,nty,ntz        ; print '(a,3i4)',     'threads per block = ',ntx,nty,ntz
      read(10,*)ibgk               ; print '(a,i1)',      'BGK order of feq  = ',ibgk
      read(10,*)ihrr               ; print '(a,i1)',      'HRR regularization= ',ihrr
      read(10,*)ivreman,smagorinsky; print '(a,i1,a,f10.4)','Vreman mixing     = ',ivreman,' Smagorinsky=',smagorinsky
      read(10,*)iforce             ; write(*,'(a,i8)',advance='no') 'iforce            = ',iforce
      select case (iforce)
      case(1)
         print '(a)','  Shan and Chen (1993)'
      case(8)
         print '(a)','  Guo (2002)'
      case(10)
         print '(a)','  Kupershtokh (2009)'
      case(12)
         print '(a)','  Khazaeli et al. 2019'
      case default
         print '(a)','  invalid forcing schemme (1,8,10,12)'
         stop
      end select
      read(10,*)nt0                ; print '(a,i8)',      'nt0               = ',nt0
      read(10,*)nt1                ; print '(a,i8)',      'nt1               = ',nt1
      if (nt1 .le. nt0) stop 'readinfile: nt1 <= nt0'
      read(10,*)iout               ; print '(a,i8)',      'iout              = ',iout
      read(10,*)irestart           ; print '(a,i8)',      'irestart          = ',irestart
      read(10,*)iprt1,iprt2,dprt   ; print '(a,3i8)',     'iprt1, iprt2, dprt= ',iprt1,iprt2,dprt
      read(10,*)lprtmin            ; print '(a,tr7,l1)',  'lprtmin           = ',lprtmin
      read(10,*)ifout              ; print '(a,i8)',      'ifout             = ',ifout
      read(10,*)ibnd               ; print '(a,i8)',      'ibnd              = ',ibnd
      read(10,*)jbnd               ; print '(a,i8)',      'jbnd              = ',jbnd
      read(10,*)kbnd               ; print '(a,i8)',      'kbnd              = ',kbnd
      read(10,*)uini,udir          ; print '(a,2(f8.3,a))','inflow (uini,udir)= ',uini,       ' [m/s]',udir,' [degrees]'
      read(10,*)rho0               ; print '(a,f8.3,a)',  'rho0 (latt dens)  = ',rho0,       ' []'
      read(10,*)rhoa               ; print '(a,f8.3,a)',  'rhoa (pres grad)  = ',rhoa,       ' []'
      read(10,'(1x,l1)')linipert   ; print '(a,tr7,l1)',  'linipert          = ',linipert
      read(10,*)inflowturbulence,turbulence_ampl,nrturb
                  ; print '(a,tr7,l1,tr2,g13.5,tr2,i5)',  'inflowturbulence  = ',inflowturbulence,turbulence_ampl,nrturb
      read(10,*)kinevisc           ; print '(a,f8.3,a)',  'Kinematic viscos  = ',kinevisc,   ' [m^2/2] '
      read(10,*)p2l%rho            ; print '(a,f8.3,a)',  'air density       = ',p2l%rho,    ' [kg/m^3]'   ! 1.225 is Air density
      read(10,*)p2l%length         ; print '(a,f8.3,a)',  'grid cell size    = ',p2l%length, ' [m]'
      read(10,*)p2l%vel            ; print '(a,f8.3,a)',  'wind velocity     = ',p2l%vel,    ' [m/s]'
      uini=uini/p2l%vel            ; print '(a,f8.3,a)',  'Non-dim uinflow   = ',uini,       ' [] Should be less that 0.2'
      read(10,'(1x,l1,1x,l1)')laveraging,laveturb ; print '(a,tr7,2l1)',  'laveraging, lavetu= ',laveraging,laveturb
      read(10,*)avestart           ; print '(a,i8)',      'avestart iteration= ',avestart
      read(10,*)avesave            ; print '(a,i8)',      'avesave iteration = ',avesave

      read(10,'(a)')ca
      if (ca /= '#-T') then
         print *,'#-T: error in infile.in'
         stop
      endif
      print *
      read(10,*)nturbines              ; print '(a,i8)',          'Num of turbines   = ',nturbines
      if (nturbines > 0) then
         allocate(ipos(nturbines), jpos(nturbines), kpos(nturbines))
         read(10,*)pitchangle          ; print '(a,f8.3,a)',      'Pitch angle       = ',pitchangle,  ' [deg]'
         read(10,*)turbrpm             ; print '(a,f8.3,a)',      'RPM for act.line  = ',turbrpm,     ' [rotations/min]'
         read(10,*)tipspeedratio       ; print '(a,f8.3,a)',      'Tipspeed ratio    = ',tipspeedratio, ' []'
         read(10,*)itiploss            ; print '(a,i8)',          'Tiploss           = ',itiploss
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
   print '(a,f12.6,a)','p2l%length = ',p2l%length,   ' [m]'
   p2l%time=p2l%length/p2l%vel
   print '(a,f12.6,a)','p2l%time   = ',p2l%time,     ' [s]'
   print '(a,f12.6,a)','p2l%vel    = ',p2l%vel,      ' [m/s]'
   print '(a,f12.6,a)','p2l%rho    = ',p2l%rho,      ' [kg/m^3]'
   p2l%visc=p2l%length**2/p2l%time
   print '(a,f12.4,a)','p2l%visc   = ',p2l%visc,     ' [m^2/s]'
   print *


!  Compute non dimensional tau from input dmensional kinematic viscosity
   print '(a,g13.6,a)',  'Kinematic visc       = ',kinevisc     ,' [m^2/s]'
   tauin = 0.5 + 3.0*kinevisc/p2l%visc
   print '(a,g13.6,a)',  'tau from kinevisc    = ',tauin        ,' [ ]'

!  Compute nondimensional kinematic viscosity used to calculate tau in fequil
   kinevisc=kinevisc/p2l%visc
   print '(a,g13.6,a)',  'Non-dim kinevisc     = ',kinevisc      ,' [ ]'

!  Compute dimensional kinematic viscosity from input non-dimensional tauin
!  tmpvisc=(1.0/3.0)*(tauin - 0.5) * p2l%visc     !(7.14)
!  print '(a,g13.6,a)',  'Kine visc from tauin = ',tmpvisc  ,' [m^2/s]'

!  Compute Reynolds number from rotor of radii lattice cells
!   reynoldsnr=(2.0*real(radii)*p2l%length)*uini*p2l%vel/(kinevisc*p2l%visc)
!   print '(a,i12,a)',    'Reynolds num         = ',nint(reynoldsnr)   ,' [ ]'

!  Compupte grid cell Reynolds number
   gridrn= p2l%length*uini*p2l%vel/(kinevisc*p2l%visc)
   print '(a,i12,a)',    'cell-Reynolds num    = ',nint(gridrn)       ,' [ ]'

! Mach number
   print '(a,f8.3,a)',   'Mach number (u/c)    = ',uini*p2l%vel/330.0 ,' [ ]'

   print *
   print '(a,g12.4)','Error terms:'
   print '(a,g12.4)','Spatial discretization errors proportional to dx^2       :', p2l%length**2
   print '(a,g12.4)','Time    discretization errors proportional to dt^2       :', p2l%time**2
   print '(a,g12.4)','Compressibility        errors proportional to dt^2/dx**2 :', p2l%time**2/p2l%length**2
   print '(a,g12.4)','BGK truncation       errors proportional to (tauin-0.5)*2:', (tauin-0.5)**2

   if (.not.runexp) stop 'runexp is false'

end subroutine
end module

