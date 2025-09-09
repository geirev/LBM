module m_turbines_forcing_kupershtokh
contains

subroutine turbines_forcing_kupershtokh(df,du,dv,dw,vel,rtmp,rho,u,v,w,dfeq1,dfeq2,ip,jp,kp,iradius,cxr,cyr,czr,nturbines,n)
!     function [U,S]=SchemeIX(A,dt,tau,f,Rho,U)
!        U= U +  0.0
!        S=( Feq(Rho,U+dt*A) - Feq(Rho,U) )./ dt
!     end
!     No update of equilibrium velocities
   use mod_dimensions
   use m_turbines_init, only : ieps
   use m_readinfile,    only : ibgk
   use m_wtime
   use mod_D3Q27setup
   use m_fequil3_block_kernel
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
   real,    intent(in) :: cxr(nl), cyr(nl), czr(nl)

   real, parameter :: inv1cs2 = 1.0/(cs2)
   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)

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
   attributes(device)  :: cxr,cyr,czr
#endif

   integer i,j,k
   integer ii
   !real !tmp1(1:100)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setting up kernel variables
#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif
   ii = 2*ieps+1

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
   !tmp1(1:3)=vel(1:3,0,63,69)
   !print '(a6,3g13.5)','vel:',!tmp1(1:3)
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute equilibrium distribution for the turbine grid points
#ifdef _CUDA
      tx=16; bx=(ii+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=4; bz=(nz+tz-1)/tz
#endif

      call fequil3_block_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(dfeq1,rtmp,vel,ii,ny,nz,nl,H2,H3,cxs,cys,czs,cs2,weights,inv1cs2,inv2cs4,inv6cs6,ibgk)
         !tmp1(1:10)=dfeq1(1:10,0,63,69)
         !print '(a6,10g13.5)','dfeq1:',!tmp1(1:10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy macro velocities plus actuator line forcing and density to vel and rtmp
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(rho, u, v, w, ip, rtmp, vel, du, dv, dw)
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
   !tmp1(1:ii)=du(-ieps:ieps,63,69)
   !print '(a6,100g13.5)','du:  ',!tmp1(1:ii)

   !tmp1(1:3)=vel(1:3,0,63,69)
   !print '(a6,3g13.5)','vel2:',!tmp1(1:3)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute equilibrium distribution for the turbine grid points with du forcing added
#ifdef _CUDA
      tx=16; bx=(ii+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=4; bz=(nz+tz-1)/tz
#endif

      call fequil3_block_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(dfeq2,rtmp,vel,ii,ny,nz,nl,H2,H3,cxs,cys,czs,cs2,weights,inv1cs2,inv2cs4,inv6cs6,ibgk)

         !tmp1(1:10)=dfeq1(1:10,0,63,69)
         !print '(a6,10g13.5)','dfeq2:',!tmp1(1:10)
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
         !tmp1(1:10)=df(1:10,0,63,69,1)
         !print '(a6,10g13.5)','df: ',!tmp1(1:10)
         !print *

end subroutine
end module
