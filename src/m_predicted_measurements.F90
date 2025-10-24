module m_predicted_measurements
type measurement
   character(len=1)c
   integer i
   integer j
   real d
end type

contains
subroutine predicted_measurements(u,v,w,it)
   use mod_dimensions
   implicit none
   real, intent(in) :: u(nx,ny,nz)
   real, intent(in) :: v(nx,ny,nz)
   real, intent(in) :: w(nx,ny,nz)
#ifdef _CUDA
   attributes(device) u,v,w
#endif
   integer, intent(in) :: it
   integer nrobs
   type(measurement), allocatable :: obs(:)
   character(len=100) fname
   character(len=5) cit
   integer iunit
   integer m,ii
   logical ex

! Reading the measurement locations and types
   inquire(file='measurement_loc.in',exist=ex)
   if (.not.ex) then
      print *,'The file measurement_loc.in does not exist'
      print *,'Are you sure lmeasurents should be true in infile,in?'
      stop
   endif

   open(newunit=iunit,file='measurement_loc.in',status='old',action='read')
      m=0
      do
        read(iunit,*,end=100,err=100)ii
        m=m+1
      enddo
      100 nrobs=m
      print *,'nrobs=',nrobs
      if (allocated(obs)) deallocate(obs); allocate(obs(nrobs))
      rewind(iunit)
      do m=1,nrobs
        read(iunit,*)obs(m)%i, obs(m)%j, obs(m)%c
      enddo
   close(iunit)

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

   write(cit,'(i5.5)')it
   fname='measurement_locations.dat'
   open(newunit=iunit,file=trim(fname), status='unknown', action='write')
      write(iunit,'(a)')'TITLE = "Measurement locations"'
      write(iunit,'(a)')'VARIABLES = "num" "i" "j" "k" "value"'
      write(iunit,'(a,i8,a)')' ZONE  F=POINT, I= ', nrobs,' J=  1, K=1'
      do m=1,nrobs
         write(iunit,'(i5,tr1,i4,i4,a,f10.5)')m,obs(m)%i,obs(m)%j,' 1 ',obs(m)%d
      enddo
   close(iunit)

end subroutine
end module


