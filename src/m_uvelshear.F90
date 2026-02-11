module m_uvelshear
contains
subroutine uvelshear(uvel_shear,uvel_time,udir_time,nt0,nt1)
! Read vertical velocity shear from file and impose on the u velocity
   use mod_dimensions
   use m_readinfile, only : udir
   implicit none
   integer, intent(in) :: nt0
   integer, intent(in) :: nt1
   real, intent(out) :: uvel_shear(nz)
   real, intent(out) :: uvel_time(nt0:nt1)
   real, intent(out) :: udir_time(nt0:nt1)
   real z(nz)
   integer k,kk
   logical ex
   real,    dimension(:),       allocatable :: uvel_h      ! temporary vertical u-velocity profile on host


! Vertical inflow velocity profile read from file
! If the uvel.dat file exists we read the vertical velocity profile into uvel(k) and normalize it
   allocate(uvel_h(nz))
   inquire(file='uvel.dat',exist=ex)
   if (ex) then
      open(10,file='uvel.dat')
         do k=1,nz
            read(10,*,err=999,end=999)kk,z(k),uvel_h(k)
         enddo
      close(10)
   else
      print '(a)','uvelshear: Did not find inputfile uvel.dat...'
      print '(a)','uvelshear: Setting constant in z velocity on inflow boundary'
      uvel_h=1.0
   endif

! Scale uvel according to uini (constant inflow velocity)
   do k=1,nz
      uvel_shear(k)=uvel_h(k)/uvel_h(nz)
   enddo
   deallocate(uvel_h)


   uvel_time(:)=1.0
   udir_time(:)=udir




   return

   999 stop 'm_uvelshear: Error reading uvel.dat'

end subroutine
end module
