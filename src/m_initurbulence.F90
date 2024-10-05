module m_initurbulence
contains
subroutine initurbulence(uu,vv,ww,rr,rho,u,v,w,inflowcor,lfirst)
   use mod_dimensions
   use m_set_random_seed2
   use m_pseudo2D
   implicit none
   real, intent(inout)  :: uu(ny,nz,0:nrturb)
   real, intent(inout)  :: vv(ny,nz,0:nrturb)
   real, intent(inout)  :: ww(ny,nz,0:nrturb)
   real, intent(inout)  :: rr(ny,nz,0:nrturb)
   real, intent(inout)  :: rho(nx,ny,nz)
   real, intent(inout)  :: u(nx,ny,nz)
   real, intent(inout)  :: v(nx,ny,nz)
   real, intent(inout)  :: w(nx,ny,nz)
   real, intent(in)     :: inflowcor
   logical, intent(in)  :: lfirst

   real cor1,cor2,dx,dy,dir
   integer(kind=4) n1,n2
   integer i,k
   integer(kind=4), save :: nxx,nyy,nzz,nt
   logical(kind=4) :: verbose=.false.

   if (lfirst) then
      call system('rm seed.dat')
      call set_random_seed2
      print '(a)','initurbulence: Computing initial turbulence field'
! Simulating smooth initial field perturbations
      cor1=20.0/sqrt(3.0)
      cor2=20.0/sqrt(3.0)
      dir=0.0
      dx=1.0
      dy=1.0
      nt=int(nrturb)
      n1=int(nx,4)
      n2=int(ny,4)
      nxx=int(nx,4)
      nyy=int(ny,4)
      nzz=int(nz,4)
      print *,'nzz=',nzz,nz
      call pseudo2D(rho,nxx,nyy,nzz,cor1,cor2,dx,dy,n1,n2,dir,verbose)
      call pseudo2D(u,nxx,nyy,nzz,cor1,cor2,dx,dy,n1,n2,dir,verbose)
      call pseudo2D(v,nxx,nyy,nzz,cor1,cor2,dx,dy,n1,n2,dir,verbose)
      call pseudo2D(w,nxx,nyy,nzz,cor1,cor2,dx,dy,n1,n2,dir,verbose)

! Imposing vertical correlation
      do k=2,nz
         rho(:,:,k)=inflowcor*rho(:,:,k-1)+sqrt(1.0-inflowcor**2)*rho(:,:,k)
         u(:,:,k)=inflowcor*u(:,:,k-1)+sqrt(1.0-inflowcor**2)*u(:,:,k)
         v(:,:,k)=inflowcor*v(:,:,k-1)+sqrt(1.0-inflowcor**2)*v(:,:,k)
         w(:,:,k)=inflowcor*w(:,:,k-1)+sqrt(1.0-inflowcor**2)*w(:,:,k)
      enddo
   endif

! Simulating a time series of inflow boundary perturbations for u
   print '(a)','initurbulence: Simulating inflow boundary turbulence'
   if (lfirst) then
      print '(a,l1)','initurbulence: lfirst=',lfirst
      uu(:,:,0)=u(1,:,:)
      vv(:,:,0)=v(1,:,:)
      ww(:,:,0)=w(1,:,:)
      rr(:,:,0)=rho(1,:,:)
   else
      print '(a,l1)','initurbulence: lfirst=',lfirst
      uu(:,:,0)=uu(:,:,nrturb)
      vv(:,:,0)=vv(:,:,nrturb)
      ww(:,:,0)=ww(:,:,nrturb)
      rr(:,:,0)=rr(:,:,nrturb)
   endif

   cor1=10.0/sqrt(3.0)
   cor2=10.0/sqrt(3.0)
   dir=0.0
   dx=1.0
   dy=1.0
   n1=int(ny,4)
   n2=int(nz,4)
   nyy=int(ny,4)
   nzz=int(nz,4)
   nt=int(nrturb)
   call pseudo2D(uu(:,:,1:nt),nyy,nzz,nt,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2D(vv(:,:,1:nt),nyy,nzz,nt,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2D(ww(:,:,1:nt),nyy,nzz,nt,cor1,cor2,dx,dy,n1,n2,dir,verbose)
   call pseudo2D(rr(:,:,1:nt),nyy,nzz,nt,cor1,cor2,dx,dy,n1,n2,dir,verbose)

! Imposing time correlations
   do i=1,nrturb
      uu(:,:,i)=inflowcor*uu(:,:,i-1)+sqrt(1.0-inflowcor**2)*uu(:,:,i)
      vv(:,:,i)=inflowcor*vv(:,:,i-1)+sqrt(1.0-inflowcor**2)*vv(:,:,i)
      ww(:,:,i)=inflowcor*ww(:,:,i-1)+sqrt(1.0-inflowcor**2)*ww(:,:,i)
      rr(:,:,i)=inflowcor*rr(:,:,i-1)+sqrt(1.0-inflowcor**2)*rr(:,:,i)
   enddo

end subroutine
end module
