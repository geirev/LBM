module m_boundarycond
contains
subroutine boundarycond(f,rho,u,v,w,fscal)
   use mod_dimensions
   use m_readinfile
   use m_bndpressure
   use m_wtime
   implicit none
   integer, parameter :: icpu=3
   real, intent(inout):: f(0:nx+1,0:ny+1,0:nz+1,nl)
   real, intent(in)   :: fscal(nl)
   real, intent(in)   :: rho(nx,ny,nz)
   real, intent(in)   :: u(nx,ny,nz)
   real, intent(in)   :: v(nx,ny,nz)
   real, intent(in)   :: w(nx,ny,nz)
   integer i,j,k,l
   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions in i-direction
   if (ibnd==0) then
! Periodic boundary conditions in i-direction
      f(0   ,:,:,:)=f(nx,:,:,:)
      f(nx+1,:,:,:)=f(1 ,:,:,:)
   elseif (ibnd==1) then
! Inflow outflow boundary conditions in i-direction
      do l=1,nl
      do k=1,nz
      do j=1,ny
         f(0,j,k,l)=fscal(l)
         f(nx+1,j,k,l)=f(nx,j,k,l)
      enddo
      enddo
      enddo
!      call bndpressure(f,rho,u,v,w)
   elseif (ibnd==2) then
! Periodic boundary conditions with pressure gradient in i-direction
      f(0   ,:,:,:)=f(nx,:,:,:)
      f(nx+1,:,:,:)=f(1 ,:,:,:)
      call bndpressure(f,rho,u,v,w)
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions in j-direction (sideways) periodic through drift.

   if (jbnd==0) then
! Periodic boundary conditions in i-direction
      f(:,0,:,:)=f(:,ny,:,:)
      f(:,ny+1,:,:)=f(:,1,:,:)
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions in k-direction (vertical)
   if (kbnd==0) then
! Periodic boundary conditions in k-direction
      f(:,:,0,:)=f(:,:,nz,:)
      f(:,:,nz+1,:)=f(: ,:,1,:)

   elseif (kbnd==1) then
! Outflow boundary on top
      f(:,:,nz+1,:)=f(:,:,nz,:)

   elseif (kbnd==2) then
! Outflow boundary on top and bottom
      f(:,:,0,:)=f(:,:,1,:)
      f(:,:,nz+1,:)=f(: ,:,nz,:)

   endif



   call cpufinish(icpu)

end subroutine
end module
