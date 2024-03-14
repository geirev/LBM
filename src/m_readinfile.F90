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

contains
subroutine readinfile
implicit none

   character(len=2) ca
   character(len=3) :: version='1.0'
   character(len=3) ver
   logical ex

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
      read(10,*)nt0                ; print '(a,tr3,i7)',   'nt0=           ',nt0
      read(10,*)nt1                ; print '(a,tr3,i7)',   'nt1=           ',nt1
      read(10,*)iout               ; print '(a,tr3,i7)',   'iout=          ',iout
      read(10,*)ifout              ; print '(a,tr3,i7)',   'ifout=         ',ifout
      read(10,*)ibnd               ; print '(a,tr3,i7)',   'ibnd=          ',ibnd
      read(10,*)jbnd               ; print '(a,tr3,i7)',   'jbnd=          ',jbnd
      read(10,*)kbnd               ; print '(a,tr3,i7)',   'kbnd=          ',kbnd
      read(10,*)uini               ; print '(a,tr2,f6.2)', 'uini=          ',uini
      read(10,*)rho0               ; print '(a,tr2,f6.2)', 'rho0=          ',rho0
      read(10,*)rhoa               ; print '(a,tr2,f6.2)', 'rhoa=          ',rhoa
      read(10,'(1x,l1)')lpseudo    ; print '(a,tr7,l1)',   'lpseudo=       ',lpseudo
      read(10,*)tau                ; print '(a,tr2,f6.2)', 'tau=           ',tau
      read(10,*)experiment         ; print '(a,tr7,a)',    'experiment=    ',trim(experiment)
   close(10)

end subroutine
end module

