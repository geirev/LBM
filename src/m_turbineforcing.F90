module m_turbineforcing
   integer, save :: iradius
   real, save :: theta=0.0
   integer, parameter :: iforce=4
contains
subroutine turbineforcing(df,rho,u,v,w,tau)
! Returns the S_i stored in df
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_fhrrscalar
   use m_actuatorline
   use m_wtime
   real, intent(out)      :: df(-ieps:ieps,ny,nz,nl,nturbines) ! forcing distributions
   real, intent(inout)    :: rho(nx,ny,nz)                     ! density
   real, intent(inout)    :: u(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: v(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: w(nx,ny,nz)                       ! velocity
   real, intent(in)       :: tau(nx,ny,nz)                     ! tau

   real cx,cy,cz,cminx,cminy,cminz,cdot,dfeq(nl)

   real                   :: force(0:ieps,ny,nz,3)             ! work array for computing the turbine force
   real                   :: ueq(-ieps:ieps,ny,nz,nturbines)   ! new turbine forced u velocity
   real                   :: veq(-ieps:ieps,ny,nz,nturbines)   ! new turbine forced v velocity
   real                   :: weq(-ieps:ieps,ny,nz,nturbines)   ! new turbine forced w velocity
   integer i,n,j,k,l,ip,jp,kp

   real, save :: dtheta=0.0
   real rps
   real A                                                      ! 0.5 for Gui and 1.0 for Kupershtokh

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

! My implementation of the actuator line method by Sørensen 2002 computing the force from all the turbines
      force=0.0

      call actuatorline(force,ny,nz,jp,kp,theta,iradius,u(ip,:,:),v(ip,:,:),w(ip,:,:),rho(ip,:,:),ieps)

! Computing the force induced velocity increments in the circular rotor plane
      if (iforce==1) then
         A=0.5
      elseif (iforce==4) then
         A=1.0
      endif

      do k=1,nz
      do j=1,ny
         if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
            do i=-ieps,ieps
               ueq(i,j,k,n)=u(ip+i,j,k)-A*force(abs(i),j,k,1)/rho(ip+i,j,k)
               veq(i,j,k,n)=v(ip+i,j,k)-A*force(abs(i),j,k,2)/rho(ip+i,j,k)
               weq(i,j,k,n)=w(ip+i,j,k)-A*force(abs(i),j,k,3)/rho(ip+i,j,k)
            enddo
         endif
      enddo
      enddo

! Computes the Gui forcing Si
      if (iforce==1) then
      ! updating the forced velocities
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  u(ip+i,j,k)=ueq(i,j,k,n)
                  v(ip+i,j,k)=veq(i,j,k,n)
                  w(ip+i,j,k)=weq(i,j,k,n)
               enddo
            endif
         enddo
         enddo

         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  do l=1,nl
                     cdot=(real(cxs(l))*ueq(i,j,k,n) + real(cys(l))*veq(i,j,k,n) + real(czs(l))*weq(i,j,k,n))/cs4

                     cminx=(real(cxs(l))-ueq(i,j,k,n))/cs2
                     cminy=(real(cys(l))-veq(i,j,k,n))/cs2
                     cminz=(real(czs(l))-weq(i,j,k,n))/cs2

                     cx=(cminx + cdot*cxs(l))*force(abs(i),j,k,1)
                     cy=(cminy + cdot*cys(l))*force(abs(i),j,k,2)
                     cz=(cminz + cdot*czs(l))*force(abs(i),j,k,3)

                     df(i,j,k,l,n)=(1.0-0.5/tau(ip+i,j,k)) * weights(l) * (cx + cy + cz)
                  enddo
               enddo
            endif
         enddo
         enddo
      endif

! Computes the new equilibrium distribution increment for the turbine planes with turbine forcing applied.
! Used for Kupershtokh forceing model.
      if (iforce==4) then
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  call fhrrscalar(dfeq,         rho(ip+i,j,k), u(ip+i,j,k), v(ip+i,j,k), w(ip+i,j,k))
                  call fhrrscalar(df(i,j,k,:,n),rho(ip+i,j,k), ueq(i,j,k,n), veq(i,j,k,n), weq(i,j,k,n))
                  df(i,j,k,:,n)=df(i,j,k,:,n)-dfeq(:)
               enddo
            else
               df(:,j,k,:,n)=0.0
            endif
         enddo
         enddo
      endif

   enddo
   call cpufinish(icpu)

end subroutine
end module

