module m_boundarycond
contains
subroutine boundarycond(f,rho,u,v,w,rr,uu,vv,ww,it,inflowvar,uvel)
   use mod_dimensions
   use m_readinfile
   use m_bndpressure
   use m_fequilscalar
   use m_wtime
   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: rho(nx,ny,nz)
   real, intent(in)   :: u(nx,ny,nz)
   real, intent(in)   :: v(nx,ny,nz)
   real, intent(in)   :: w(nx,ny,nz)
   real, intent(in)   :: uu(ny,nz,0:nrturb)
   real, intent(in)   :: vv(ny,nz,0:nrturb)
   real, intent(in)   :: ww(ny,nz,0:nrturb)
   real, intent(in)   :: rr(ny,nz,0:nrturb)
   real, intent(in)   :: inflowvar
   real, intent(in)   :: uvel(nz)
   integer, intent(in):: it
   real :: rtmp(ny,nz)
   real :: utmp(ny,nz)
   real :: vtmp(ny,nz)
   real :: wtmp(ny,nz)
   integer j,k,lit,ja,ka
   integer, parameter :: icpu=8

   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions in i-direction
   if (ibnd==0) then
! Periodic boundary conditions in i-direction
      f(:,0   ,:,:)=f(:,nx,:,:)
      f(:,nx+1,:,:)=f(:,1 ,:,:)

   elseif (ibnd==1) then
! Inflow outflow boundary conditions in i-direction
      lit=mod(it,nrturb)
      if (lit==0) lit=nrturb
      do k=1,nz
      do j=1,ny
         utmp(j,k)=uvel(k)+inflowvar*uu(j,k,lit)
         vtmp(j,k)=inflowvar*vv(j,k,lit)
         wtmp(j,k)=inflowvar*ww(j,k,lit)
         rtmp(j,k)=rho(1,j,k)+inflowvar*rr(j,k,lit)
      enddo
      enddo

      do k=0,nz+1
      ka=min(max(k,1),nz)
      do j=0,ny+1
         ja=min(max(j,1),ny)
         f(1:nl,0,j,k)=fequilscalar(rtmp(ja,ka),utmp(ja,ka),vtmp(ja,ka),wtmp(ja,ka))
      enddo
      enddo


      do k=0,nz+1
      do j=0,ny+1
         f(:,nx+1,j,k)=f(:,nx,j,k)
      enddo
      enddo

   elseif (ibnd==2) then
! Periodic boundary conditions with pressure gradient in i-direction
      f(:,0   ,:,:)=f(:,nx,:,:)
      f(:,nx+1,:,:)=f(:,1 ,:,:)
      call bndpressure(f,rho,u,v,w)
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions in j-direction (sideways) periodic through drift.

   if (jbnd==0) then
! Periodic boundary conditions in i-direction
      f(:,:,0,:)   =f(:,:,ny,:)
      f(:,:,ny+1,:)=f(:,:,1,:)
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions in k-direction (vertical)
   if (kbnd==0) then
! Periodic boundary conditions in k-direction
      f(:,:,:,0)   =f(:,:,:,nz)
      f(:,:,:,nz+1)=f(:,: ,:,1)

   elseif (kbnd==1) then
! Outflow boundary on top
      f(:,:,:,nz+1)=f(:,:,:,nz)

   elseif (kbnd==2) then
! Outflow boundary on top and bottom
      f(:,:,:,0)   =f(:,:,:,1)
      f(:,:,:,nz+1)=f(:,:,:,nz)
   endif



   call cpufinish(icpu)

end subroutine
end module
