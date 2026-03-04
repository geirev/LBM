module m_inflow
contains
subroutine inflow(uvel_shear,uvel_time,udir_time,nt0,nt1)
! Read vertical velocity shear from file and impose on the u velocity
   use mod_dimensions
   use m_readinfile, only : uini,udir,p2l
   implicit none
   integer, intent(in) :: nt0
   integer, intent(in) :: nt1
   real, intent(out) :: uvel_shear(nz)
   real, intent(out) :: uvel_time(nt0:nt1)
   real, intent(out) :: udir_time(nt0:nt1)
   real z(nz)
   integer i,k,kk
   logical ex
   real tmp,t
   integer nrtdata
   real, allocatable :: tdata(:),uvel_tdata(:),udir_tdata(:)
   real, parameter   :: pi=3.1415927410125732
   real,    dimension(:),       allocatable :: uvel_h      ! temporary vertical u-velocity profile on host

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vertical inflow velocity profile read from file
! If the uvel.dat file exists we read the vertical velocity profile into uvel(k), normalize it
! and store scaled uvel in uvel_shear
   allocate(uvel_h(nz))
   uvel_h=1.0
   inquire(file='uvel_shear.dat',exist=ex)
   if (ex) then
      print '(a)','inflow: Reading inflow vertical profile from uvel_shear.dat'
      open(10,file='uvel_shear.dat')
         do k=1,nz
            read(10,*,err=999,end=999)kk,z(k),uvel_h(k)
         enddo
      close(10)
   endif

   do k=1,nz
      uvel_shear(k)=uvel_h(k)/uvel_h(nz)
   enddo
   deallocate(uvel_h)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Time variability of inflow-velocity and inflow-direction
! First reading a time series of times, velocities and directions from uvel_time.dat
! Next interpolate to continuous time series in uvel_time(nt0:nt1) and udir_time(nt0:nt1)
   uvel_time(:)=uini
   udir_time(:)=udir
   inquire(file='uvel_time.dat',exist=ex)
   if (ex) then
      print '(a)','inflow: Reading inflow time variability from uvel_time.dat'
      open(10,file='uvel_time.dat')
         nrtdata = 0
         do
            read(10,*,err=998,end=100) tmp
            nrtdata = nrtdata + 1
         end do
   100   rewind(10)

         allocate(tdata(nrtdata),uvel_tdata(nrtdata),udir_tdata(nrtdata))
         do k=1,nrtdata
            read(10,*,err=998,end=998) tdata(k),uvel_tdata(k),udir_tdata(k)
            uvel_tdata(k)=uvel_tdata(k)/p2l%vel
         enddo
      close(10)

      do i=nt0,nt1
         t = real(i-1)*p2l%time
         if (t <= tdata(1)) then
            uvel_time(i) = uvel_tdata(1)
            udir_time(i) = udir_tdata(1)
         else if (t >= tdata(nrtdata)) then
            uvel_time(i) = uvel_tdata(nrtdata)
            udir_time(i) = udir_tdata(nrtdata)
         else
            do k=1,nrtdata-1
               if (t >= tdata(k) .and. t <= tdata(k+1)) then
                  tmp = (t - tdata(k)) / (tdata(k+1) - tdata(k))
                  uvel_time(i) = uvel_tdata(k) + tmp*(uvel_tdata(k+1) - uvel_tdata(k))
                  udir_time(i) = udir_tdata(k) + tmp*(udir_tdata(k+1) - udir_tdata(k))
                  exit
               endif
            enddo
         endif
      enddo
   endif

!   do k=nt0,nt1
!      udir_time(k)=0.0 + 20.0*sin(real(k)*pi*2.0/real(nt1-nt0))
!   enddo


   return

   998 stop 'm_inflow: Error reading uvel_time.dat'
   999 stop 'm_inflow: Error reading uvel_shear.dat'

end subroutine
end module
