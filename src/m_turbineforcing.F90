module m_turbineforcing
   integer, save :: iradius
   real, save :: theta=0.0

   real, allocatable, public   :: turbine_df(:,:,:,:,:)      ! Turbine forcing used in applyturbines
   real, allocatable, private  :: dfeq1(:,:,:,:)
   real, allocatable, private  :: dfeq2(:,:,:,:)
   real, allocatable, private  :: force(:,:,:,:)     ! work array for computing the turbine force
   real, allocatable, private  :: du(:,:,:)          ! turbine forced u velocity
   real, allocatable, private  :: dv(:,:,:)          ! turbine forced v velocity
   real, allocatable, private  :: dw(:,:,:)          ! turbine forced w velocity
   real, allocatable, private  :: vel(:,:,:,:)
   real, allocatable, private  :: rtmp(:,:,:)
   real, allocatable, private  :: cx(:), cy(:), cz(:)
   real, allocatable, private  :: slice_u(:,:)
   real, allocatable, private  :: slice_v(:,:)
   real, allocatable, private  :: slice_w(:,:)
   real, allocatable, private  :: slice_r(:,:)
#ifdef _CUDA
   attributes(device) :: turbine_df
   attributes(device) :: dfeq1
   attributes(device) :: dfeq2
   attributes(device) :: force
   attributes(device) :: du
   attributes(device) :: dv
   attributes(device) :: dw
   attributes(device) :: vel
   attributes(device) :: rtmp
   attributes(device) :: cx, cy, cz
   attributes(device) :: slice_u
   attributes(device) :: slice_v
   attributes(device) :: slice_w
   attributes(device) :: slice_r
#endif

contains
subroutine init_turbines
   use mod_dimensions
   use m_readinfile, only : nturbines
   use mod_D3Q27setup
   implicit none
   integer l
   if (.not. allocated(turbine_df))  allocate(turbine_df(1:nl,-ieps:ieps,1:ny,1:nz,nturbines))
   if (.not. allocated(dfeq1))       allocate(dfeq1(nl,-ieps:ieps,ny,nz)            )
   if (.not. allocated(dfeq2))       allocate(dfeq2(nl,-ieps:ieps,ny,nz)            )
   if (.not. allocated(force))       allocate(force(0:ieps,ny,nz,3))
   if (.not. allocated(du))          allocate(du(-ieps:ieps,ny,nz) )
   if (.not. allocated(dv))          allocate(dv(-ieps:ieps,ny,nz) )
   if (.not. allocated(dw))          allocate(dw(-ieps:ieps,ny,nz) )
   if (.not. allocated(vel))         allocate(vel(3,-ieps:ieps,ny,nz) )
   if (.not. allocated(rtmp))        allocate(rtmp(-ieps:ieps,ny,nz) )
   if (.not. allocated(cx))          allocate(cx(nl))
   if (.not. allocated(cy))          allocate(cy(nl))
   if (.not. allocated(cz))          allocate(cz(nl))
   if (.not. allocated(slice_u))     allocate(slice_u(ny,nz))
   if (.not. allocated(slice_v))     allocate(slice_v(ny,nz))
   if (.not. allocated(slice_w))     allocate(slice_w(ny,nz))
   if (.not. allocated(slice_r))     allocate(slice_r(ny,nz))

   do l=1,nl
     cx(l)=real(cxs_h(l))
     cy(l)=real(cys_h(l))
     cz(l)=real(czs_h(l))
   enddo
end subroutine

subroutine turbineforcing(rho,u,v,w,it,nt1)
! Returns the S_i stored in turbine_df and possibly updated velocities
   use mod_dimensions, only : nx,ny,nz,ieps
   use mod_D3Q27setup
   use m_readinfile,   only : turbrpm,p2l,ipos,jpos,kpos,nturbines
   use m_fequilscalar
   use m_actuatorline
   use m_turbineforcing_kupershtokh
#ifdef _CUDA
   use cudafor
#endif
   use m_wtime
   real, intent(inout) :: rho(nx,ny,nz)                     ! density
   real, intent(inout) :: u(nx,ny,nz)                       ! velocity
   real, intent(inout) :: v(nx,ny,nz)                       ! velocity
   real, intent(inout) :: w(nx,ny,nz)                       ! velocity
#ifdef _CUDA
   attributes(device)  :: rho
   attributes(device)  :: u
   attributes(device)  :: v
   attributes(device)  :: w
#endif
   !!!!real cminx,cminy,cminz,cdotuu

   real :: force_h(0:ieps,ny,nz,3)        ! work array for computing the turbine force

   integer :: it,nt1
   integer i,n,j,k,ip,jp,kp
   real :: u_h(ny, nz)
   real :: v_h(ny, nz)
   real :: w_h(ny, nz)
   real :: r_h(ny, nz)

   real, save :: dtheta=0.0
   real rps

   real, parameter :: pi=3.1415927410125732
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0
   integer, parameter :: icpu=2
   call cpustart()

! Rotations per timestep
! Starting with turbrpm (rotations per minute)
   rps=turbrpm/60.0    ! rotations per second

! rotations per nondim time step
   dtheta=rps*p2l%time*pi2


   theta=theta+dtheta

   turbine_df(:,:,:,:,:)=0.0
   du(:,:,:)=0.0
   dv(:,:,:)=0.0
   dw(:,:,:)=0.0


   do n=1,nturbines
      ip=ipos(n)
      jp=jpos(n)
      kp=kpos(n)

! My implementation of the actuator line method by Sørensen 2002 computing the force from all the turbines
      force_h=0.0

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k = 1, nz
        do j = 1, ny
          slice_u(j,k) = u(ip,j,k)
          slice_v(j,k) = v(ip,j,k)
          slice_w(j,k) = w(ip,j,k)
          slice_r(j,k) = rho(ip,j,k)
        end do
      end do

      u_h = slice_u
      v_h = slice_v
      w_h = slice_w
      r_h = slice_r

      call actuatorline(force_h,ny,nz,jp,kp,theta,iradius,u_h,v_h,w_h,r_h,ieps)

      force(0:ieps,1:ny,1:nz,1:3)=force_h(0:ieps,1:ny,1:nz,1:3)


! Computing the force induced velocity increments in the circular 3D rotor plane (F/rho)
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(tubine_du, tubine_dv, tubine_dw, tubine_force, rho, ip, jp, kp, iradius)
#endif
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
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

! (10) Kupershtokh 2009
      call  turbineforcing_kupershtokh(turbine_df, du, dv, dw, vel, rtmp, rho, u, v, w, dfeq1, dfeq2,&
                                       ip, jp, kp, iradius, cx, cy, cz, nturbines, n, it, nt1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! (8) Guo 2002
!     function [U,S]=SchemeVII(A,dt,tau,f,Rho,U)
!        U= U +  dt*A./2
!        S= Rho.*W’.*(Cdot(A)- sum(A.*(U+du),1))./ Cs2 *(1-1/(2*tau))
!        S=S+Rho.*W’.*(Cdot(A) .* Cdot(U+du) )./ Cs2^2*(1-1/(2*tau))
!     end
!!        if (iforce==8) then
!!  #ifdef _CUDA
!!  !$cuf kernel do(2) <<<*,*>>>
!!  #endif
!!           do k=1,nz
!!           do j=1,ny
!!              if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
!!                 do i=-ieps,ieps
!!                    ! updating the forced equilibrium velocities ueq=u + F/(2rho)
!!                    u(ip+i,j,k)=u(ip+i,j,k)+0.5*du(i,j,k)
!!                    v(ip+i,j,k)=v(ip+i,j,k)+0.5*dv(i,j,k)
!!                    w(ip+i,j,k)=w(ip+i,j,k)+0.5*dw(i,j,k)
!!  
!!                    ! Computing the S_i term returned in df
!!  !                  cdota(:)=real(cxs(:))*du(i,j,k) + real(cys(:))*dv(i,j,k) + real(czs(:))*dw(i,j,k)
!!  !                  cdotu(:)=real(cxs(:))*u(ip+i,j,k) + real(cys(:))*v(ip+i,j,k) + real(czs(:))*w(ip+i,j,k)
!!  !                  udota   =u(ip+i,j,k)*du(i,j,k) + v(ip+i,j,k)*dv(i,j,k) + w(ip+i,j,k)*dw(i,j,k)
!!  !
!!  !                  df(i,j,k,:,n)= rho(ip+i,j,k) * weights(:) * ( (cdota(:) - udota)/cs2   + (cdota(:) * cdotu(:) )/cs4 )
!!  
!!                    do l=1,nl
!!                       cdotuu=(real(cxs(l))*u(ip+i,j,k) + real(cys(l))*v(ip+i,j,k) + real(czs(l))*w(ip+i,j,k))/cs4
!!                       cminx=(real(cxs(l))-u(ip+i,j,k))/cs2
!!                       cminy=(real(cys(l))-v(ip+i,j,k))/cs2
!!                       cminz=(real(czs(l))-w(ip+i,j,k))/cs2
!!  
!!                       cx=(cminx + cdotuu*real(cxs(l)))*du(i,j,k)
!!                       cy=(cminy + cdotuu*real(cys(l)))*dv(i,j,k)
!!                       cz=(cminz + cdotuu*real(czs(l)))*dw(i,j,k)
!!  
!!                       df(l,i,j,k,n)= rho(ip+i,j,k) * weights(l) * (cx + cy + cz)
!!                    enddo
!!                 enddo
!!              endif
!!           enddo
!!           enddo
!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!      elseif (iforce==10) then

!!        elseif (iforce==12) then       
!!  !(12) Khazaeli et al. 2019
!!  !     function [U,S]=SchemeXI(A,dt,tau,f,Rho,U)
!!  !        S=(1-1/(2*tau))*( Feq(Rho,U+dt*A) - Feq(Rho,U) )./ dt
!!  !        U= U +  dt*A./2
!!  !     end
!!  #ifdef _CUDA
!!  !$cuf kernel do(2) <<<*,*>>>
!!  #endif
!!           do k=1,nz
!!           do j=1,ny
!!              if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
!!                 do i=-ieps,ieps
!!  !                  dfeq(:)      =fequilscalar(rho(ip+i,j,k), u(ip+i,j,k),           v(ip+i,j,k),           w(ip+i,j,k))
!!  !                  df(:,i,j,k,n)=fequilscalar(rho(ip+i,j,k), u(ip+i,j,k)+du(i,j,k), v(ip+i,j,k)+dv(i,j,k), w(ip+i,j,k)+dw(i,j,k))
!!                    df(:,i,j,k,n)= (df(:,i,j,k,n)-dfeq(:))
!!  
!!                    u(ip+i,j,k)=u(ip+i,j,k)+0.5*du(i,j,k)
!!                    v(ip+i,j,k)=v(ip+i,j,k)+0.5*dv(i,j,k)
!!                    w(ip+i,j,k)=w(ip+i,j,k)+0.5*dw(i,j,k)
!!                 enddo
!!              endif
!!           enddo
!!           enddo
!!  
!!  
!!        else
!!           stop 'turbine forcing uses invalid forcing model iforce'
!!        endif

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
