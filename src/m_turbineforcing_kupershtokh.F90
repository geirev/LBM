module m_turbineforcing_kupershtokh
contains

subroutine turbineforcing_kupershtokh(df,du,dv,dw,vel,rtmp,&
                                      rho,u,v,w,dfeq1,dfeq2,ip,jp,kp,iradius,cx,cy,cz,nturbines,n,it,nt1)
!     function [U,S]=SchemeIX(A,dt,tau,f,Rho,U)
!        U= U +  0.0
!        S=( Feq(Rho,U+dt*A) - Feq(Rho,U) )./ dt
!     end
!     No update of equilibrium velocities
   use mod_dimensions
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
   integer :: it,nt1,i1,i2

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

   integer i,j,k,ii

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setting up kernel variables
#ifdef _CUDA
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
