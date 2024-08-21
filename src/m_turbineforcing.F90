module m_turbineforcing
contains
subroutine turbineforcing(df,feq,rho,u,v,w)
   use mod_dimensions
   use m_readinfile
   use m_feqscalar
   use m_actuatorline
   use m_wtime
   real, intent(out)      :: df(ny,nz,nl,nturbines) ! forcing distributions
   real, intent(in)       :: feq(0:nx+1,0:ny+1,0:nz+1,nl) ! equilibrium distribution
   real, intent(in)       :: rho(nx,ny,nz)          ! density
   real, intent(in)       :: u(nx,ny,nz)            ! velocity
   real, intent(in)       :: v(nx,ny,nz)            ! velocity
   real, intent(in)       :: w(nx,ny,nz)            ! velocity

   real                   :: force(ny,nz,3)         ! work array for computing the turbine force
   integer n,j,k,ip,jp,kp

   integer, save :: iradius
   real, save :: theta=0.0
   real, save :: dtheta=0.0
   real rps,rpts

   real, parameter :: pi=3.1415926535
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0
   integer, parameter :: icpu=5
   call cpustart()


! Rotations per timestep
! Starting with turbrpm (rotations per minute)
   rps=turbrpm/60.0    ! rotations per second
   rpts=rps*p2l%time   ! rotations per time step of p2l%time
   dtheta=rpts*pi2


   theta=theta+dtheta
   do n=1,nturbines
      ip=ipos(n)
      jp=jpos(n)
      kp=kpos(n)
      force=0.0

! My implementation of the actuator line method by Sørensen 2002 computing the force from all the turbines
      call actuatorline(force(1:ny,1:nz,1:3),ny,nz,jp,kp,theta,turbrpm,iradius,u(ip,:,:),v(ip,:,:),w(ip,:,:))

! Computes the new equilibrium distribution for the turbine planes with turbine forcing applied.
      df(:,:,:,n)=0.0
      do k=1,nz
      do j=1,ny
         if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
            call feqscalar(df(j,k,:,n),rho(ip,j,k),u(ip,j,k)-force(j,k,1),v(ip,j,k)-force(j,k,2),w(ip,j,k)-force(j,k,3))
            !call feqscalar(df(j,k,:,n),rho(ip,j,k),u(ip,j,k)-force(j,k,1),v(ip,j,k),w(ip,j,k))
         else
            df(j,k,:,n)=feq(ip,j,k,:)
         endif
      enddo
      enddo

! The returned forcing distribution (will be zero outside the turbine planes).
      df(1:ny,1:nz,1:nl,n)=df(1:ny,1:nz,1:nl,n)-feq(ip,1:ny,1:nz,1:nl)
   enddo
   call cpufinish(icpu)

end subroutine
end module

