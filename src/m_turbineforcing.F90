module m_turbineforcing
   integer, save :: iradius
   real, save :: theta=0.0
contains
subroutine turbineforcing(df,rho,u,v,w,tau)
! Returns the S_i stored in df and possibly updated velocities
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

   real cx,cy,cz,cminx,cminy,cminz,cdotu,dfeq(nl)

   real                   :: force(0:ieps,ny,nz,3)         ! work array for computing the turbine force
   real                   :: du(-ieps:ieps,ny,nz)          ! turbine forced u velocity
   real                   :: dv(-ieps:ieps,ny,nz)          ! turbine forced v velocity
   real                   :: dw(-ieps:ieps,ny,nz)          ! turbine forced w velocity
   integer i,n,j,k,l,ip,jp,kp

   integer, parameter :: iforce=1
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

! Used for Guo (2002) model
      if (iforce==1) then
      ! updating the forced equilibrium velocities ueq=u + F/(2rho)
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  u(ip+i,j,k)=u(ip+i,j,k)+0.5*du(i,j,k)
                  v(ip+i,j,k)=v(ip+i,j,k)+0.5*dv(i,j,k)
                  w(ip+i,j,k)=w(ip+i,j,k)+0.5*dw(i,j,k)
               enddo
            endif
         enddo
         enddo

      ! Computing the S_i term returned in df
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  do l=1,nl
                     cdotu=(real(cxs(l))*u(ip+i,j,k) + real(cys(l))*v(ip+i,j,k) + real(czs(l))*w(ip+i,j,k))/cs4

                     cminx=(real(cxs(l))-u(ip+i,j,k))/cs2
                     cminy=(real(cys(l))-v(ip+i,j,k))/cs2
                     cminz=(real(czs(l))-w(ip+i,j,k))/cs2

                     cx=(cminx + cdotu*real(cxs(l)))*force(abs(i),j,k,1)
                     cy=(cminy + cdotu*real(cys(l)))*force(abs(i),j,k,2)
                     cz=(cminz + cdotu*real(czs(l)))*force(abs(i),j,k,3)

                     df(i,j,k,l,n)=(1.0-0.5/tau(ip+i,j,k)) * weights(l) * (cx + cy + cz)
                  enddo
               enddo
            endif
         enddo
         enddo
      endif

! Used for Kupershtokh forcing model.
      if (iforce==4) then
      ! No update of equilibrium velocities

      ! Computing the S_i term returned in df
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                  dfeq(:)      =fhrrscalar(rho(ip+i,j,k), u(ip+i,j,k),           v(ip+i,j,k),           w(ip+i,j,k))
                  df(i,j,k,:,n)=fhrrscalar(rho(ip+i,j,k), u(ip+i,j,k)+du(i,j,k), v(ip+i,j,k)+dv(i,j,k), w(ip+i,j,k)+dw(i,j,k))
                  df(i,j,k,:,n)=df(i,j,k,:,n)-dfeq(:)

                  !df(i,j,k,:,n)=fhrrscalar(rho(ip+i,j,k), u(ip+i,j,k)+du(i,j,k), v(ip+1,j,k)+dv(i,j,k), w(ip+1,j,k)+dw(i,j,k))&
                  !             -fhrrscalar(rho(ip+i,j,k), u(ip+i,j,k)          , v(ip+i,j,k)          , w(ip+i,j,k))
               enddo
            endif
         enddo
         enddo
      endif

   enddo
   call cpufinish(icpu)

end subroutine
end module

!-    for it=1:MaxIter
!- ! Calculate Macroscopic Variables
!-       [Rho,U]=CalculateMacroscopic(f);
!-
!- ! Calculate Macroscopic Force (Acceleration)
!-       [A]=CalculateAcceleration(Rho,U);
!-
!- ! Calculate Microscopic Force (S)
!- ! This function name “ForceModel” should be replaced with the proper name of the desired force scheme.
!-       [du,S,Uph]=ForceModel(A,dt,tau,f,Rho,U);
!-
!- ! Apply Microscopic Force to velocity
!-       Ueq=U+du; ! calculate the equilibrium velocity
!-
!- ! Perform the collision
!-       [f]=Collision(tau,f,Rho,Ueq);
!-
!- ! Apply Microscopic Force to f
!-       f=f+dt*S; !adding force
!- ! Perform the streaming
!-
!-       f=Streaming(f);
!- ! Apply the boundary conditions
!-       f=Boundaries(f,Rho,Ueq);
!-    end
!-
!-
!- ! A Helper Function
!- ! Dot Lattice speed with variable X
!-    function [res]=Cdot(X)
!-       sizes=size(x);
!-       sizes(1)=1;
!-       res=zeros(sizes);
!-       for id=1:size(X,1)
!-          res = res +C(id,:)’.* X(id,:);
!-       end
!-    end
!-
!- ! Shan and Chen 1993
!-    function [du,S,Uph]=SchemeI(A,dt,tau,f,Rho,U)
!-       du=dt*tau*A;
!-       Uph=U+A.*dt./2;
!-       S=zeros(size(f));
!-    end
!-
!- ! Luo 1997
!-    function [du,S, Uph]=SchemeII(A,dt,tau,f,Rho,U)
!-       du= zeros(size(U));
!-       Uph=U;
!-       S=W’ .* Cdot(Rho.*A) ./ Cs2;
!-    end
!-
!- !  He et al., 1998
!-    function [du,S, Uph]=SchemeIII(A,dt,tau,f,Rho,U)
!-       du= zeros(size(U));
!-       Uph=U;
!-       feq=Feq(Rho,U);
!-       S=feq .* (Cdot(A)- sum(A.*(U),1))./ Cs2;
!-    end
!-
!- ! Martys 1998 and Luo 1998
!-    function [du,S, Uph]=SchemeIV(A,dt,tau,f,Rho,U)
!-       du= zeros(size(U));
!-       Uph=U;
!-       S= Rho.*W’.*(Cdot(A)- sum(A.*U,1))./ Cs2;
!-       S=S+Rho.*W’.*(Cdot(A) .* Cdot(U) )./ Cs2^2;
!-    end
!-
!- ! Buick and Greated 2000
!-    function [du,S, Uph]=SchemeV(A,dt,tau,f,Rho,U)
!-       du= dt*A./2;
!-       Uph=U+du;
!-       S=(1-1/(2*tau))*W’.*Cdot(Rho.*A)./ Cs2;
!-    end
!-
!- ! Ladd and Verberg 2001
!-    function [du,S, Uph]=SchemeVIa(A,dt,tau,f,Rho,U)
!-       du= zeros(size(U));
!-       Uph=U;
!-       S= Rho.*W’.*(Cdot(A)- sum(A.*U,1))./ Cs2;
!-       S=S+Rho.*W’.*(Cdot(A) .* Cdot(U) )./ Cs2^2*(1-1/(2*tau));
!-       S=S+Rho.*W’.* Cdot(A) ./ Cs2 /(2*tau);
!-    end
!-
!- ! Ladd and Verberg 2001
!-    function [du,S, Uph]=SchemeVIb(A,dt,tau,f,Rho,U)
!-       du= zeros(size(U));
!-       Uph=U+dt*A./2;
!-       S= Rho.*W’.*(Cdot(A)- sum(A.*U,1))./ Cs2;
!-       S=S+Rho.*W’.*(Cdot(A) .* Cdot(U) )./ Cs2^2;
!-    end
!-
!- ! Guo 2002
!-    function [du,S, Uph]=SchemeVII(A,dt,tau,f,Rho,U)
!-       du= dt*A./2;
!-       Uph=U+du;
!-       S= Rho.*W’.*(Cdot(A)- sum(A.*(U+du),1))./ Cs2 *(1-1/(2*tau));
!-       S=S+Rho.*W’.*(Cdot(A) .* Cdot(U+du) )./ Cs2^2*(1-1/(2*tau));
!-    end
!-
!- ! Wagner 2006
!-    function [du,S, Uph]=SchemeVIII(A,dt,tau,f,Rho,U)
!-       du= zeros(size(U));
!-       Uph=U+dt*A./2;
!-       S= Rho.*W’.*(Cdot(A)- sum(A.*U,1))./ Cs2 ;
!-       S=S+Rho.*W’.*(Cdot(A) .* Cdot(U) )./ Cs2^2;
!-       S=S+Rho.*W’.*(Cdot(A) .* Cdot(A) )./ Cs2^2*(1-1/(4*tau))*dt/2;
!-       S=S-Rho.*W’.*( sum(A.*A,1))./ Cs2 *(1-1/(4*tau))*dt/2;
!-    end
!-
!- ! Kupershtokh 2009
!-    function [du,S, Uph]=SchemeIX(A,dt,tau,f,Rho,U)
!-       du= zeros(size(U));
!-       Uph=U+dt*A./2;
!-       S=( Feq(Rho,U+dt*A) - Feq(Rho,U) )./ dt;
!-    end
!-
!- ! Guo et al. 2011
!-    function [du,S, Uph]=SchemeX(A,dt,tau,f,Rho,U)
!-       du= +dt*A./2;
!-       Uph=U+du;
!-       feq=Feq(Rho,U+du);
!-       S=feq .* (Cdot(A)- sum(A.*(U+du),1))./ Cs2*(1-1/(2*tau));
!-    end
!-
!- ! Khazaeli et al. 2019
!-    function [du,S, Uph]=SchemeXI(A,dt,tau,f,Rho,U)
!-       du= dt*A./2;
!-       Uph=U+dt*A./2;
!-       S=( Feq(Rho,U+dt*A) - Feq(Rho,U) )./ dt*(1-1/(2*tau));
!-    end
