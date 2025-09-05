!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

module m_turbines_forcing_guo
contains

subroutine turbines_forcing_guo(df,du,dv,dw,vel,rtmp,&
                                      rho,u,v,w,dfeq1,dfeq2,ip,jp,kp,iradius,cx,cy,cz,nturbines,n,it,nt1)
! (8) Guo 2002
!     function [U,S]=SchemeVII(A,dt,tau,f,Rho,U)
!        U= U +  dt*A./2
!        S= Rho.*W’.*(Cdot(A)- sum(A.*(U+du),1))./ Cs2 *(1-1/(2*tau))
!        S=S+Rho.*W’.*(Cdot(A) .* Cdot(U+du) )./ Cs2^2*(1-1/(2*tau))
!     end
   use mod_dimensions
   use m_turbines_init, only : ieps
   use m_wtime
   use mod_D3Q27setup
   use m_fequilscal
   use m_fequilscalar
#ifdef _CUDA
   use cudafor
#endif
   implicit none

! Computing the S_i term returned in df
   integer, intent(in) :: nturbines
   real, intent(inout) :: df(nl,-ieps:ieps,ny,nz,nturbines)
   real, intent(in)    :: du(-ieps:ieps,ny,nz)
   real, intent(in)    :: dv(-ieps:ieps,ny,nz)
   real, intent(in)    :: dw(-ieps:ieps,ny,nz)

   real, intent(in)    :: rho(nx,ny,nz)
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)

   real, intent(inout) :: vel(3,-ieps:ieps,ny,nz)
   real, intent(inout) :: rtmp(-ieps:ieps,ny,nz)

   real, intent(inout) :: dfeq1(nl,-ieps:ieps,ny,nz)
   real, intent(inout) :: dfeq2(nl,-ieps:ieps,ny,nz)

   integer, intent(in) :: ip, jp, kp, iradius, n
   real,    intent(in) :: cx(nl), cy(nl), cz(nl)
   integer :: it,nt1

#ifdef _CUDA
   attributes(device)  :: df
   attributes(device)  :: du
   attributes(device)  :: dv
   attributes(device)  :: dw
   attributes(device)  :: rho
   attributes(device)  :: u
   attributes(device)  :: v
   attributes(device)  :: w
   attributes(device)  :: vel
   attributes(device)  :: rtmp
   attributes(device)  :: dfeq1
   attributes(device)  :: dfeq2
   attributes(device)  :: cx,cy,cz
#endif

   integer i,j,k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setting up kernel variables
#ifdef _CUDA
   integer ii
   type(dim3) :: grid,tblock

   ii = 2*ieps+1

   tBlock%x = 1   ! for i
   tBlock%y = 8   ! for j
   tBlock%z = 8   ! for k

   grid%x = (ii + tBlock%x - 1) / tBlock%x
   grid%y = (ny + tBlock%y - 1) / tBlock%y
   grid%z = (nz + tBlock%z - 1) / tBlock%z
#endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy macro velocities and density to rtmp and vel

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(rho, u, v, w, ip, rtmp, vel)
#endif
   do k=1,nz
   do j=1,ny
   do i=-ieps,ieps
      rtmp(i,j,k) =rho(ip+i,j,k)
      vel(1,i,j,k)=u(ip+i,j,k)
      vel(2,i,j,k)=v(ip+i,j,k)
      vel(3,i,j,k)=w(ip+i,j,k)
   enddo
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute equilibrium distribution for the turbine grid points

#ifdef _CUDA
   call fequilscal<<<grid, tBlock>>>(dfeq1, rtmp, vel, weights, cx, cy, cz, H2, H3, ii)
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(jp,kp, ieps, iradius, dfeq1, rtmp, vel, weights, cx, cy, cz, H2, H3)
   do k=1,nz
   do j=1,ny
      if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
         do i=-ieps,ieps
            dfeq1(:,i,j,k)=fequilscalar(rtmp(i,j,k), vel(1,i,j,k), vel(2,i,j,k), vel(3,i,j,k), weights, cx, cy, cz, H2, H3)
         enddo
      endif
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy macro velocities plus actuator line forcing and density to vel and rtmp

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(rho, u, v, w, ip, rtmp, vel, du, dv, dw, ip)
#endif
   do k=1,nz
   do j=1,ny
   do i=-ieps,ieps
      rtmp(i,j,k) =rho(ip+i,j,k)
      vel(1,i,j,k)=u(ip+i,j,k)+du(i,j,k)
      vel(2,i,j,k)=v(ip+i,j,k)+dv(i,j,k)
      vel(3,i,j,k)=w(ip+i,j,k)+dw(i,j,k)
   enddo
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute equilibrium distribution for the turbine grid points with du forcing added

#ifdef _CUDA
   call fequilscal<<<grid, tBlock>>>(dfeq2, rtmp, vel, weights, cx, cy, cz, H2, H3, ii)
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(jp, kp, ieps, dfeq2, rtmp, vel, weights, cx, cy, cz, H2, H3 )
   do k=1,nz
   do j=1,ny
      if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
         do i=-ieps,ieps
            dfeq2(:,i,j,k)=fequilscalar(rtmp(i,j,k), vel(1,i,j,k), vel(2,i,j,k), vel(3,i,j,k), weights, cx, cy, cz, H2, H3)
         enddo
      endif
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute final turbine_df to be used in applyturbines

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(df, dfeq1, dfeq2, jp, kp, iradius )
#endif
         do k=1,nz
         do j=1,ny
            if ( ((j-jp)**2 + (k-kp)**2 ) <  (iradius+5)**2) then
               do i=-ieps,ieps
                   df(:,i,j,k,n)=dfeq2(:,i,j,k)-dfeq1(:,i,j,k)
               enddo
            endif
         enddo
         enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif


end subroutine
end module
