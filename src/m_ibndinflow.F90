module m_ibndinflow
contains
subroutine ibndinflow(feqscal,u)
   use mod_dimensions
   use m_feqscalar
   use m_readinfile
   implicit none
   real, intent(out) :: feqscal(0:nz+1,nl)
   real, intent(out) :: u(nx,ny,nz)

   real z(nz)
   real uvel(nz)
   integer i,j,k,kk
   logical ex

   inquire(file='uvel.dat',exist=ex)
   if (.not.ex) then
      print '(a)','ibndinflow: Did not find inputfile uvel.dat...'
      stop
   endif
   open(10,file='uvel.dat')
      do k=1,nz
         read(10,*,err=999,end=999)kk,z(k),uvel(k)
      enddo
   close(10)

   do k=1,nz
      uvel(k)=uini*uvel(k)/uvel(nz)
   enddo

   do k=1,nz
      call feqscalar(feqscal(k,:),rho0,uvel(k),0.0,0.0)
   enddo

   feqscal(0,:)=feqscal(1,:)
   feqscal(nz+1,:)=feqscal(nz,:)

   do j=1,ny
   do i=1,nx
      u(i,j,1:nz)=uvel(1:nz)
   enddo
   enddo
   return
   999 stop 'ibndinflow: Error reading uvel.dat'

end subroutine
end module
