module m_turbineforcing
   integer, save :: iradius
   real, save :: theta=0.0
contains
subroutine turbineforcing(df,feq,rho,u,v,w)
   use mod_dimensions
   use m_readinfile
   use m_fhrrscalar
   use m_actuatorline
   use m_wtime
   real, intent(out)      :: df(-ieps:ieps,ny,nz,nl,nturbines)     ! forcing distributions
   real, intent(in)       :: feq(0:nx+1,0:ny+1,0:nz+1,nl) ! equilibrium distribution
   real, intent(in)       :: rho(nx,ny,nz)          ! density
   real, intent(in)       :: u(nx,ny,nz)            ! velocity
   real, intent(in)       :: v(nx,ny,nz)            ! velocity
   real, intent(in)       :: w(nx,ny,nz)            ! velocity

   real                   :: force(0:ieps,ny,nz,3)         ! work array for computing the turbine force
   integer i,n,j,k,ip,jp,kp

   real, save :: dtheta=0.0
   real rps

   real, parameter :: pi=3.1415927410125732
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0
   integer, parameter :: icpu=5
   call cpustart()


! Rotations per timestep
! Starting with turbrpm (rotations per minute)
   rps=turbrpm/60.0    ! rotations per second

! rotations per nondim time step
   dtheta=rps*p2l%time*pi2


   theta=theta+dtheta
   df(:,:,:,:,:)=0.0
   do n=1,nturbines
      ip=ipos(n)
      jp=jpos(n)
      kp=kpos(n)
      force=0.0

! My implementation of the actuator line method by Sørensen 2002 computing the force from all the turbines
      call actuatorline(force(0:ieps,1:ny,1:nz,1:3),ny,nz,jp,kp,theta,iradius,u(ip,:,:),v(ip,:,:),w(ip,:,:),rho(ip,:,:),ieps)

! Computes the new equilibrium distribution for the turbine planes with turbine forcing applied.
      do k=1,nz
      do j=1,ny
         if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
            do i=-ieps,ieps
               call fhrrscalar(df(i,j,k,:,n),rho(ip+i,j,k),&
                                            u(ip+i,j,k)-force(abs(i),j,k,1)/rho(ip+i,j,k),&
                                            v(ip+i,j,k)-force(abs(i),j,k,2)/rho(ip+i,j,k),&
                                            w(ip+i,j,k)-force(abs(i),j,k,3)/rho(ip+i,j,k))
            enddo
         else
            do i=-ieps,ieps
               df(i,j,k,:,n)=feq(ip+i,j,k,:)
            enddo
         endif
      enddo
      enddo

! The returned forcing distribution (will be zero outside the turbine planes).
      df(-ieps:ieps,1:ny,1:nz,1:nl,n)=df(-ieps:ieps,1:ny,1:nz,1:nl,n)-feq(ip-ieps:ip+ieps,1:ny,1:nz,1:nl)
   enddo
   call cpufinish(icpu)

end subroutine
end module

