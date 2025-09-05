module m_predicted_measurements
logical :: lmeasurements=.true.
type measurement
   character(len=1)c
   integer i
   integer j
   real d
end type

contains
subroutine predicted_measurements(u,v,w,it)
   use mod_dimensions
   real, intent(in) :: u(nx,ny,nz)
   real, intent(in) :: v(nx,ny,nz)
   real, intent(in) :: w(nx,ny,nz)
#ifdef _CUDA
   attributes(device) u,v,w
#endif
   integer, intent(in) :: it
   integer, parameter :: nrobs=6
   type(measurement) obs(nrobs)
   character(len=100) fname
   character(len=5) cit
   integer iunit

   obs(1)%i=70 ; obs(1)%j=28; obs(1)%c='u'
   obs(2)%i=70 ; obs(2)%j=28; obs(2)%c='v'
   obs(3)%i=86 ; obs(3)%j=85; obs(3)%c='u'
   obs(4)%i=86 ; obs(4)%j=85; obs(4)%c='v'
   obs(5)%i=102; obs(5)%j=62; obs(5)%c='u'
   obs(6)%i=102; obs(6)%j=62; obs(6)%c='v'

   do m=1,nrobs
      if (obs(m)%c == 'u') then
         obs(m)%d=u(obs(m)%i,obs(m)%j,1)
      elseif (obs(m)%c == 'v') then
         obs(m)%d=v(obs(m)%i,obs(m)%j,1)
      else
         print *,'invalid variable identifier'
      endif
   enddo

   write(cit,'(i5.5)')it
   fname='measurements_'//cit//'.dat'
   open(newunit=iunit,file=trim(fname), status='unknown', action='write')
      do m=1,nrobs
         write(iunit,'(i5,tr1,a1,i4,i4,f10.5)')m,obs(m)%c,obs(m)%i,obs(m)%j,obs(m)%d
      enddo
   close(iunit)

end subroutine
end module
