module m_turbineforcing
contains
subroutine turbineforcing(df,feq,rho,u,v,w)
   use mod_dimensions
   use m_readinfile
   use m_feqscalar
   use m_actuatorline
   real, intent(out)      :: df(ny,nz,nl,nturbines) ! forcing distributions
   real, intent(in)       :: feq(0:nx+1,0:ny+1,0:nz+1,nl) ! equilibrium distribution
   real, intent(in)       :: rho(nx,ny,nz)          ! density
   real, intent(in)       :: u(nx,ny,nz)            ! velocity
   real, intent(in)       :: v(nx,ny,nz)            ! velocity
   real, intent(in)       :: w(nx,ny,nz)            ! velocity

   real                   :: force(ny,nz)           ! work array for computing the turbine force
   integer n,j,k,ip,jp,kp

   real, save :: theta=0.0
   real, save :: dtheta=0.0
   real rps,rpts

   real, parameter :: pi=3.1415926535
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0


! Rotations per timestep
   rps=turbrpm/60.0
   rpts=rps*p2l%time
   dtheta=rpts*pi2


   theta=theta+dtheta

      do n=1,nturbines
         ip=ipos(n)
         jp=jpos(n)
         kp=kpos(n)
         force=0.0
         call actuatorline(force(1:ny,1:nz),ny,nz,jp,kp,real(radii),theta)
         force=1.0-force

         df(:,:,:,n)=0.0
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  radii**2) then
               call feqscalar(df(j,k,:,n),rho(ip,j,k),force(j,k)*u(ip,j,k),force(j,k)*v(ip,j,k),force(j,k)*w(ip,j,k))
            else
               df(j,k,:,n)=feq(ip,j,k,:)
            endif
         enddo
         enddo
         df(1:ny,1:nz,1:nl,n)=df(1:ny,1:nz,1:nl,n)-feq(ip,1:ny,1:nz,1:nl)
      enddo

!      do n=1,nturbines
!         ip=ipos(n)
!         jp=jpos(n)
!         kp=kpos(n)
!         df(:,:,:,n)=0.0
!         do k=1,nz
!         do j=1,ny
!            if ( ((j-jp)**2 + (k-kp)**2 ) <  radii**2) then
!               call feqscalar(df(j,k,:,n),rho(ip,j,k),turbblock*u(ip,j,k),turbblock*v(ip,j,k),turbblock*w(ip,j,k))
!            else
!               df(j,k,:,n)=feq(ip,j,k,:)
!            endif
!         enddo
!         enddo
!         df(1:ny,1:nz,1:nl,n)=df(1:ny,1:nz,1:nl,n)-feq(ip,1:ny,1:nz,1:nl)
!      enddo

end subroutine
end module

