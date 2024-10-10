module m_channeldiag
contains
subroutine channeldiag(f,rho,u,lblanking)
   use mod_dimensions
   use m_readinfile, only : rho0,rhoa,tauin
   use m_density
   use m_velocity
   use mod_D3Q27setup
   implicit none
   real, intent(inout)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)    :: u(nx,ny,nz)
   real, intent(inout)    :: rho(nx,ny,nz)
   logical, intent(inout) :: lblanking(nx,ny,nz)
   real width,mu,grad,dw
   integer j

   rho=density(f,lblanking)
   u= velocity(f,rho,cxs,lblanking)
   width=real(ny-2)
   mu=cs2*(tauin-0.5)
   grad=cs2*2.0*rhoa/real(nx)
   dw=width/real(ny-1)
   open(10,file='channeluvel.dat')
   do j=1,ny
      write(10,'(i4,2f13.5)')j,u(nx/2,j,2),&
                 (1.0/(2.0*rho0*mu))*grad*((width/2.0)**2-(dw*(real(j)-1.0)-width/2.0)**2)
   enddo
   close(10)
end subroutine
end module
