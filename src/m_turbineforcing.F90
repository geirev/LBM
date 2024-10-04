module m_turbineforcing
   integer, save :: iradius
   real, save :: theta=0.0
contains
subroutine turbineforcing(df,rho,u,v,w,tau)
! Returns the S_i stored in df and possibly updated velocities
   use mod_dimensions, only : nx,ny,nz,ieps
   use mod_D3Q27setup
   use m_readinfile,  only : turbrpm,p2l,ipos,jpos,kpos,nturbines,iforce
   use m_fequilscalar
   use m_actuatorline
   use m_wtime
   real, intent(out)      :: df(-ieps:ieps,ny,nz,nl,nturbines) ! forcing distributions
   real, intent(inout)    :: rho(nx,ny,nz)                     ! density
   real, intent(inout)    :: u(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: v(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: w(nx,ny,nz)                       ! velocity
   real, intent(in)       :: tau(nx,ny,nz)                     ! tau

   real cx,cy,cz,cminx,cminy,cminz,cdotuu,dfeq(nl)
   real cdotA(nl),udotA,cdotu(nl)

   real                   :: force(0:ieps,ny,nz,3)         ! work array for computing the turbine force
   real                   :: du(-ieps:ieps,ny,nz)          ! turbine forced u velocity
   real                   :: dv(-ieps:ieps,ny,nz)          ! turbine forced v velocity
   real                   :: dw(-ieps:ieps,ny,nz)          ! turbine forced w velocity
   integer i,n,j,k,l,ip,jp,kp

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

! My implementation of the actuator line method by Sørensen 2002 computing the force from all the turbines
      force=0.0

      call actuatorline(force,ny,nz,jp,kp,theta,iradius,u(ip,:,:),v(ip,:,:),w(ip,:,:),rho(ip,:,:),ieps)

! Computing the force induced velocity increments in the circular 3D rotor plane (F/rho)
      do k=1,nz
      do j=1,ny
         if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
            do i=-ieps,ieps
               du(i,j,k)=-force(abs(i),j,k,1)/rho(ip+i,j,k)
               dv(i,j,k)=-force(abs(i),j,k,2)/rho(ip+i,j,k)
               dw(i,j,k)=-force(abs(i),j,k,3)/rho(ip+i,j,k)
            enddo
         endif
      enddo
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! (1) Shan and Chen 1993
!     function [U,S]=SchemeI(A,dt,tau,f,Rho,U)
!        U= U + dt*tau*A
!        S=0.0
!     end
      if (iforce==1) then
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  ! updating the forced equilibrium velocities
                  u(ip+i,j,k)=u(ip+i,j,k)+tau(ip+i,j,k)*du(i,j,k)
                  v(ip+i,j,k)=v(ip+i,j,k)+tau(ip+i,j,k)*dv(i,j,k)
                  w(ip+i,j,k)=w(ip+i,j,k)+tau(ip+i,j,k)*dw(i,j,k)

                  ! No update of S_i term

               enddo
            endif
         enddo
         enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! (8) Guo 2002
!     function [U,S]=SchemeVII(A,dt,tau,f,Rho,U)
!        U= U +  dt*A./2
!        S= Rho.*W’.*(Cdot(A)- sum(A.*(U+du),1))./ Cs2 *(1-1/(2*tau))
!        S=S+Rho.*W’.*(Cdot(A) .* Cdot(U+du) )./ Cs2^2*(1-1/(2*tau))
!     end
      elseif (iforce==8) then
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  ! updating the forced equilibrium velocities ueq=u + F/(2rho)
                  u(ip+i,j,k)=u(ip+i,j,k)+0.5*du(i,j,k)
                  v(ip+i,j,k)=v(ip+i,j,k)+0.5*dv(i,j,k)
                  w(ip+i,j,k)=w(ip+i,j,k)+0.5*dw(i,j,k)

                  ! Computing the S_i term returned in df
!                  cdota(:)=real(cxs(:))*du(i,j,k) + real(cys(:))*dv(i,j,k) + real(czs(:))*dw(i,j,k)
!                  cdotu(:)=real(cxs(:))*u(ip+i,j,k) + real(cys(:))*v(ip+i,j,k) + real(czs(:))*w(ip+i,j,k)
!                  udota   =u(ip+i,j,k)*du(i,j,k) + v(ip+i,j,k)*dv(i,j,k) + w(ip+i,j,k)*dw(i,j,k)
!
!                  df(i,j,k,:,n)=(1.0-0.5/tau(ip+i,j,k)) * rho(ip+i,j,k) * weights(:) * &
!                        ( (cdota(:) - udota)/cs2   + (cdota(:) * cdotu(:) )/cs4 )

                  do l=1,nl
                     cdotuu=(real(cxs(l))*u(ip+i,j,k) + real(cys(l))*v(ip+i,j,k) + real(czs(l))*w(ip+i,j,k))/cs4
                     cminx=(real(cxs(l))-u(ip+i,j,k))/cs2
                     cminy=(real(cys(l))-v(ip+i,j,k))/cs2
                     cminz=(real(czs(l))-w(ip+i,j,k))/cs2

                     cx=(cminx + cdotuu*real(cxs(l)))*du(i,j,k)
                     cy=(cminy + cdotuu*real(cys(l)))*dv(i,j,k)
                     cz=(cminz + cdotuu*real(czs(l)))*dw(i,j,k)

                     df(i,j,k,l,n)=(1.0-0.5/tau(ip+i,j,k)) * rho(ip+i,j,k) * weights(l) * (cx + cy + cz)
                  enddo
               enddo
            endif
         enddo
         enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! (10) Kupershtokh 2009
!     function [U,S]=SchemeIX(A,dt,tau,f,Rho,U)
!        U= U +  0.0
!        S=( Feq(Rho,U+dt*A) - Feq(Rho,U) )./ dt
!     end
      elseif (iforce==10) then
      ! No update of equilibrium velocities

      ! Computing the S_i term returned in df
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  dfeq(:)      =fequilscalar(rho(ip+i,j,k), u(ip+i,j,k),           v(ip+i,j,k),           w(ip+i,j,k))
                  df(i,j,k,:,n)=fequilscalar(rho(ip+i,j,k), u(ip+i,j,k)+du(i,j,k), v(ip+i,j,k)+dv(i,j,k), w(ip+i,j,k)+dw(i,j,k))
                  df(i,j,k,:,n)=df(i,j,k,:,n)-dfeq(:)
               enddo
            endif
         enddo
         enddo

!(12) Khazaeli et al. 2019
!     function [U,S]=SchemeXI(A,dt,tau,f,Rho,U)
!        S=(1-1/(2*tau))*( Feq(Rho,U+dt*A) - Feq(Rho,U) )./ dt
!        U= U +  dt*A./2
!     end
      elseif (iforce==12) then
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  dfeq(:)      =fequilscalar(rho(ip+i,j,k), u(ip+i,j,k),           v(ip+i,j,k),           w(ip+i,j,k))
                  df(i,j,k,:,n)=fequilscalar(rho(ip+i,j,k), u(ip+i,j,k)+du(i,j,k), v(ip+i,j,k)+dv(i,j,k), w(ip+i,j,k)+dw(i,j,k))
                  df(i,j,k,:,n)= (df(i,j,k,:,n)-dfeq(:))   ! Moved the tau part to apply turbines since tau is recomputed in feq

                  u(ip+i,j,k)=u(ip+i,j,k)+0.5*du(i,j,k)
                  v(ip+i,j,k)=v(ip+i,j,k)+0.5*dv(i,j,k)
                  w(ip+i,j,k)=w(ip+i,j,k)+0.5*dw(i,j,k)
               enddo
            endif
         enddo
         enddo


      else
         stop 'turbine forcing uses invalid forcing model iforce'
      endif

   enddo
   call cpufinish(icpu)

end subroutine
end module

!    for it=1:MaxIter
!       [Rho,U]= CalculateMacroscopic(f)
!       [A]    = CalculateAcceleration(Rho,U)
!       [Ueq,S]= ForceModel(A,dt,tau,f,Rho,U)
!       [f]    = Collision(tau,f,Rho,Ueq)
!       f      = f+dt*S; !adding force
!       f      = Streaming(f)
!       f      = Boundaries(f,Rho,Ueq)
!    end
!
!
!   Dot Lattice speed with variable X
!    function [res(1:nl)]=Cdot(X)
!       res(:) = res(:) + cxs(:)*X(1) + cys(:)*X(2)  + czs(:)*X(3)
!    end
!
!   (1) Shan and Chen 1993
!    function [U,S]=SchemeI(A,dt,tau,f,Rho,U)
!       du=dt*tau*A
!       U= U + du
!       S=0
!    end
!
!   (2) Luo 1997
!    function [U,S]=SchemeII(A,dt,tau,f,Rho,U)
!       du=0
!       U= U + du
!       S=W’ * Cdot(Rho*A) / Cs2
!    end
!
!   (3)  He et al., 1998
!    function [U,S]=SchemeIII(A,dt,tau,f,Rho,U)
!       du=0
!       U= U + du
!       S=Feq(Rho,U) * (Cdot(A)- sum(A*(U),1)) / Cs2
!    end
!
!   (4) Martys 1998 and Luo 1998
!    function [U,S]=SchemeIV(A,dt,tau,f,Rho,U)
!       du=0
!       U= U + du
!       S= Rho*W’*(Cdot(A)- sum(A*U,1))/ Cs2
!       S=S+Rho*W’*(Cdot(A) * Cdot(U) )/ Cs2^2
!    end
!
!   (5) Buick and Greated 2000
!    function [U,S]=SchemeV(A,dt,tau,f,Rho,U)
!       du=dt*A/2
!       U= U + du
!       S=(1-1/(2*tau))*W’*Cdot(Rho*A)/ Cs2
!    end
!
!   (6) Ladd and Verberg 2001
!    function [U,S]=SchemeVIa(A,dt,tau,f,Rho,U)
!       du=0
!       U= U + du
!       S= Rho*W’*(Cdot(A)- sum(A*U,1))/ Cs2
!       S=S+Rho*W’*(Cdot(A) * Cdot(U) )/ Cs2^2*(1-1/(2*tau))
!       S=S+Rho*W’* Cdot(A) / Cs2 /(2*tau)
!    end
!
!   (7) Ladd and Verberg 2001
!    function [U,S]=SchemeVIb(A,dt,tau,f,Rho,U)
!       du=0
!       U= U + du
!       S= Rho*W’*(Cdot(A)- sum(A*U,1))/ Cs2
!       S=S+Rho*W’*(Cdot(A) * Cdot(U) )/ Cs2^2
!    end
!
!   (8) Guo 2002
!    function [U,S]=SchemeVII(A,dt,tau,f,Rho,U)
!       du=dt*A/2
!       U= U + du
!       S= (Cdot(A)- sum(A*(U+du),1))/ Cs2
!       S=S + (Cdot(A) * Cdot(U+du) )/ Cs2^2
!       S=(1-1/(2*tau))*Rho*W’*S
!    end
!
!   (9) Wagner 2006
!    function [U,S]=SchemeVIII(A,dt,tau,f,Rho,U)
!       du=0
!       U= U + du
!       S= Rho*W’*(Cdot(A)- sum(A*U,1))/ Cs2
!       S=S+Rho*W’*(Cdot(A) * Cdot(U) )/ Cs2^2
!       S=S+Rho*W’*(Cdot(A) * Cdot(A) )/ Cs2^2*(1-1/(4*tau))*dt/2
!       S=S-Rho*W’*( sum(A*A,1))/ Cs2 *(1-1/(4*tau))*dt/2
!    end
!
!   (10) Kupershtokh 2009
!    function [U,S]=SchemeIX(A,dt,tau,f,Rho,U)
!       du=0
!       U= U + du
!       S=( Feq(Rho,U+dt*A) - Feq(Rho,U) )/ dt
!    end
!
!   (11) Guo et al. 2011
!    function [U,S]=SchemeX(A,dt,tau,f,Rho,U)
!       du=dt*A/2
!       U= U + du
!       S=Feq(Rho,U+du) * (Cdot(A)- sum(A*(U+du),1))/ Cs2*(1-1/(2*tau))
!    end
!
!   (12) Khazaeli et al. 2019
!    function [U,S]=SchemeXI(A,dt,tau,f,Rho,U)
!       du=dt*A/2
!       S=( Feq(Rho,U+dt*A) - Feq(Rho,U) )/ dt*(1-1/(2*tau))
!       U= U + du
!    end
