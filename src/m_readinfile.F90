module m_readinfile
! Simulation parameters
   integer  nt0      ! First timestep
   integer  nt1      ! Last timestep
   integer  iout     ! number of steps between outputs 0, 0+iout, ...
   real     rho0     ! Average density
   real     tau      ! Collision timescale 0.6
   logical  dbg      ! Print diagnostics if tru
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
      read(10,*)nt0                ; print '(a,tr3,i5)',   'nt0=           ',nt0
      read(10,*)nt1                ; print '(a,tr3,i5)',   'nt1=           ',nt1
      read(10,*)iout               ; print '(a,tr3,i5)',   'iout=          ',iout
      read(10,*)rho0               ; print '(a,tr2,f6.2)', 'rho0=          ',rho0
      read(10,*)tau                ; print '(a,tr2,f6.2)', 'tau=           ',tau
      read(10,'(1x,l1)')dbg        ; print '(a,tr7,l1)',   'dbg=           ',dbg
      read(10,*)experiment         ; print '(a,tr7,a)',    'experiment=    ',trim(experiment)
   close(10)

end subroutine
end module

