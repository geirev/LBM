module m_readinfile
! Simulation parameters
   integer  nt0            ! First timestep
   integer  nt1            ! Last timestep
   integer  iout           ! number of steps between outputs 0, 0+iout, ...
   integer  iprt1          ! Output every dprt time steps of it <= iprt
   integer  iprt2          ! Output every dprt time steps of it <= iprt
   integer  dprt           ! delta high frequency output 
   logical  ltesting       ! Print minimalistice plt file if true (no derived variables)
   integer  irestart       ! number of steps between restart files
   integer  itecout        ! format of tecplot solution files (0 full files, 2 only solution variables)
   integer  ibnd           ! Type of bondary condition in i direction (ibnd=0 periodic, 1 inflow/outflow, 12
   integer  jbnd           ! Type of bondary condition in i direction
   integer  kbnd           ! Type of bondary condition in k direction
   logical  inflowturbulence          ! Add smooth pseudo-random peturbations for turbulent inflow
   integer  nrturb         ! Number of precomputed batches of inflow turbulence for u,v,w, and rho
   real     turbulence_ampl! strength of inflow turbulence
   real     uini           ! Initial absolute velocity
   real     udir           ! Initial velocity direction in degrees
   real     rho0           ! Average density
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
   logical :: ldump        ! Dumping diagnostic files to disk
   logical :: lmeasurements! Used in data assimilation experiments

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
   real, allocatable ::  yaw(:),tilt(:)          ! Turbine yaw and tilt
   integer  ihrr           ! Option (1) for regularized R(fneq) scheme
   integer  ibgk           ! Option (2,3) for second or third order BGK f^eq expansion
   integer  ivreman        ! Option (1) for subgridscale mixing using Vreman
   integer  iablvisc       ! Atmospheric boundary layer mixing (0-none, 1-mechanical layer, 2-buoyancy scheme)
   real     ablheight      ! Height of atmospheric boundary layer
   integer  istable        ! stability of abl (1=stable, 0=neutral, and -1=unstable)
   real smagorinsky        ! smagorinsky constant (0.15) used in subgridscale mixing
   logical ltiming         ! true for timing of kernels. Should be false to avoid syncs

contains
subroutine readinfile()
   use m_mkinfile
   implicit none

   character(len=1) ver
   logical ex
   real gridrn
   integer n

   inquire(file='main.F90',exist=ex)
   if (ex) stop 'You are executing boltzmann in the src/ catalog'

! reading input data
   inquire(file='infile.in',exist=ex)
   if (.not.ex) then
      print '(a)','Did not find inputfile infile.in'
      print '(a)','Generating new template infile.in.....'
      call mkinfile()
      print '(a)','Please edit and relaunch boltzmann'
      stop
   endif
   print '(a)','--------------------------------------------------------------------------------'
   open(10,file='infile.in')
      read(10,'(a)',err=100)ver

      read(10,*,err=100)ltiming            ; print '(a,tr7,l1)',  'ltiming           = ',ltiming
      read(10,*,err=100)ltesting           ; print '(a,tr7,l1)',  'ltesting          = ',ltesting
      read(10,*,err=100)ldump              ; print '(a,tr7,l1)',  'ldump             = ',ldump
      read(10,*,err=100)ntx,nty,ntz        ; print '(a,3i4)',     'threads per block = ',ntx,nty,ntz
      read(10,*,err=100)lmeasurements      ; print '(a,tr7,l1)',  'lmeasurements     = ',lmeasurements

      read(10,'(a)',err=100)ver

      read(10,*,err=100)experiment         ; print '(a,a)',       'experiment        = ',trim(experiment)
      read(10,*,err=100)ibgk               ; print '(a,i1)',      'BGK order of feq  = ',ibgk
      read(10,*,err=100)ihrr               ; print '(a,i1)',      'HRR regularization= ',ihrr
      read(10,*,err=100)ivreman,smagorinsky; print '(a,i1,a,f10.4)','Vreman mixing     = ',ivreman,' Smagorinsky=',smagorinsky

      read(10,'(a)',err=100)ver

      read(10,*,err=100)iablvisc           ; print '(a,i1)',      'ABL: iablvisc     = ',iablvisc
      read(10,*,err=100)ablheight          ; print '(a,f10.4)',   'ABL: ablheight    = ',ablheight
      read(10,*,err=100)istable            ; print '(a,i3)',      'ABL: stability    = ',istable

      read(10,'(a)',err=100)ver

      read(10,*,err=100)nt0                ; print '(a,i8)',      'nt0               = ',nt0
      read(10,*,err=100)nt1                ; print '(a,i8)',      'nt1               = ',nt1
      if (nt1 .le. nt0) stop 'readinfile: nt1 <= nt0'
      read(10,*,err=100)iout               ; print '(a,i8)',      'iout              = ',iout
      read(10,*,err=100)irestart           ; print '(a,i8)',      'irestart          = ',irestart
      read(10,*,err=100)iprt1,iprt2,dprt   ; print '(a,3i8)',     'iprt1, iprt2, dprt= ',iprt1,iprt2,dprt
      read(10,*,err=100)itecout            ; print '(a,i8)',      'itecout           = ',itecout

      read(10,'(a)',err=100)ver

      read(10,*,err=100)ibnd               ; print '(a,i8)',      'ibnd              = ',ibnd
      read(10,*,err=100)jbnd               ; print '(a,i8)',      'jbnd              = ',jbnd
      read(10,*,err=100)kbnd               ; print '(a,i8)',      'kbnd              = ',kbnd

      read(10,'(a)',err=100)ver

      read(10,*,err=100)uini,udir          ; print '(a,2(f8.3,a))','inflow (uini,udir)= ',uini,       ' [m/s]',udir,' [degrees]'
      read(10,*,err=100)inflowturbulence,turbulence_ampl,nrturb
                  ; print '(a,tr7,l1,tr2,g13.5,tr2,i5)',  'inflowturbulence  = ',inflowturbulence,turbulence_ampl,nrturb

      read(10,'(a)',err=100)ver

      read(10,*,err=100)kinevisc           ; print '(a,f8.3,a)',  'Kinematic viscos  = ',kinevisc,   ' [m^2/2] '
      read(10,*,err=100)p2l%rho            ; print '(a,f8.3,a)',  'air density       = ',p2l%rho,    ' [kg/m^3]'   ! 1.225 is Air density
      read(10,*,err=100)p2l%length         ; print '(a,f8.3,a)',  'grid cell size    = ',p2l%length, ' [m]'
      read(10,*,err=100)p2l%vel            ; print '(a,f8.3,a)',  'wind velocity     = ',p2l%vel,    ' [m/s]'
      uini=uini/p2l%vel            ; print '(a,f8.3,a)',  'Non-dim uinflow   = ',uini,       ' [] Should be less that 0.2'

      read(10,'(a)',err=100)ver

      read(10,'(1x,l1,1x,l1)',err=100)laveraging,laveturb ; print '(a,tr7,2l1)',  'laveraging, lavetu= ',laveraging,laveturb
      read(10,*,err=100)avestart           ; print '(a,i8)',      'avestart iteration= ',avestart
      read(10,*,err=100)avesave            ; print '(a,i8)',      'avesave iteration = ',avesave

      read(10,'(a)',err=100)ver

      read(10,*,err=100)nturbines              ; print '(a,i8)',          'Num of turbines   = ',nturbines
      if (nturbines > 0) then
         allocate(ipos(nturbines), jpos(nturbines), kpos(nturbines), yaw(n), tilt(n))
         read(10,*,err=100)pitchangle          ; print '(a,f8.3,a)',      'Pitch angle       = ',pitchangle,  ' [deg]'
         read(10,*,err=100)turbrpm             ; print '(a,f8.3,a)',      'RPM for act.line  = ',turbrpm,     ' [rotations/min]'
         read(10,*,err=100)tipspeedratio       ; print '(a,f8.3,a)',      'Tipspeed ratio    = ',tipspeedratio, ' []'
         read(10,*,err=100)itiploss            ; print '(a,i8)',          'Tiploss           = ',itiploss
!         do n=1,nturbines
!            read(10,'(a)',err=100)ver
!            read(10,*,err=100)ipos(n)
!            read(10,*,err=100)jpos(n)
!            read(10,*,err=100)kpos(n)
!            print '(a,i4,a,3i4)', '(ijk)-pos for turbine  = ',n,' : ',ipos(n),jpos(n),kpos(n)
!         enddo
         do n=1,nturbines
            read(10,'(a)',err=100)ver
            read(10,*,err=100)ipos(n),jpos(n),kpos(n),yaw(n),tilt(n)
            print '(a,i4,a,3i4,2f10.2)', '(ijk)-pos for turbine  = ',n,' : ',ipos(n),jpos(n),kpos(n),yaw(n),tilt(n)
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


!  Compute non dimensional tau from input dimensional kinematic viscosity
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

   rho0=1.0
   return

   100 close(10)
   call system('mv -i infile.in infile_backup.in')
   call mkinfile()
   print *,'infile.in problem; infile.in moved to infile_backup.in and generated new template infile.in'
   stop

end subroutine
end module

