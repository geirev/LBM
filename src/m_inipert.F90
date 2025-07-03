module m_inipert
contains
subroutine inipert(rho,u,v,w,uvel)
   use mod_dimensions
   use m_set_random_seed2
   use m_pseudo2D
   use m_readinfile, only : rho0,linipert,udir
   implicit none
   real, intent(inout)  :: rho(nx,ny,nz)
   real, intent(inout)  :: u(nx,ny,nz)
   real, intent(inout)  :: v(nx,ny,nz)
   real, intent(inout)  :: w(nx,ny,nz)
   real, intent(in)     :: uvel(nz)

!   real :: vertcor=0.95
   real :: stddev=0.00001  ! value running stable in uniform flow for 2000 timesteps
   real, parameter :: pi=3.1415927410125732

!   real cor1,cor2,dx,dy,dir
!   integer(kind=4) n1,n2
   integer k
!   integer(kind=4) :: nxx,nyy,nzz
!   logical(kind=4) :: verbose=.false.

   u=0.0
   v=0.0
   w=0.0
   rho=rho0

   if (linipert) then
! Simulating smooth initial field perturbations
!      cor1=20.0/sqrt(3.0)
!      cor2=20.0/sqrt(3.0)
!      dir=0.0
!      dx=1.0
!      dy=1.0
!      n1=int(nx,4)
!      n2=int(ny,4)
!      nxx=int(nx,4)
!      nyy=int(ny,4)
!      nzz=int(nz,4)
!      call pseudo2D(rho,nxx,nyy,nzz,cor1,cor2,dx,dy,n1,n2,dir,verbose)
!      call pseudo2D(u,nxx,nyy,nzz,cor1,cor2,dx,dy,n1,n2,dir,verbose)
!      call pseudo2D(v,nxx,nyy,nzz,cor1,cor2,dx,dy,n1,n2,dir,verbose)
!      call pseudo2D(w,nxx,nyy,nzz,cor1,cor2,dx,dy,n1,n2,dir,verbose)

! Imposing vertical correlation
!      do k=2,nz
!         rho(:,:,k)=vertcor*rho(:,:,k-1)+sqrt(1.0-vertcor**2)*rho(:,:,k)
!         u(:,:,k)=vertcor*u(:,:,k-1)+sqrt(1.0-vertcor**2)*u(:,:,k)
!         v(:,:,k)=vertcor*v(:,:,k-1)+sqrt(1.0-vertcor**2)*v(:,:,k)
!         w(:,:,k)=vertcor*w(:,:,k-1)+sqrt(1.0-vertcor**2)*w(:,:,k)
!      enddo

      call random_number(rho)
      rho=rho0 + stddev*rho
      do k=1,nz
         u(:,:,k)=uvel(k)*cos(udir*pi/180.0)
         v(:,:,k)=uvel(k)*sin(udir*pi/180.0)
      enddo
      w=0.0
   endif

end subroutine
end module
