module m_uvelshear
contains
subroutine uvelshear(uvel)
! Read vertical velocity shear from file and impose on the u velocity
   use mod_dimensions
   use m_readinfile, only : uini
   implicit none
   real, intent(out) :: uvel(nz)
   real z(nz)
   integer k,kk
   logical ex

! If the uvel.dat file exists we read the normalized vertical velocity profile into uvel(k)
   inquire(file='uvel.dat',exist=ex)
   if (ex) then
      open(10,file='uvel.dat')
         do k=1,nz
            read(10,*,err=999,end=999)kk,z(k),uvel(k)
         enddo
      close(10)
   else
      print '(a)','ibndinflow: Did not find inputfile uvel.dat...'
      print '(a)','ibndinflow: Setting constant in z velocity on inflow boundary'
      uvel=1.0
   endif

! Scale uvel according to uini (constant inflow velocity)
   do k=1,nz
      uvel(k)=uini*uvel(k)/uvel(nz)
   enddo

   return
   999 stop 'ibndinflow: Error reading uvel.dat'

end subroutine
end module
