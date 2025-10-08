module m_inflow_turbulence_forcing
! kupershtokh forcing is used for the turbulence.
! Its also applied on a thin slice of one grid cell thickness.
contains
subroutine inflow_turbulence_forcing(rho,u,v,w,ampl,it,nrturb)
   use mod_dimensions
   use m_inflow_turbulence_init
   use m_fequil3_block_kernel
   use m_readinfile,    only : ibgk
#ifdef _CUDA
   use cudafor
#endif
   use mod_D3Q27setup
   use m_wtime

   integer, intent(in)    :: nrturb
   real, intent(inout)    :: rho(nx,ny,nz)                     ! density
   real, intent(inout)    :: u(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: v(nx,ny,nz)                       ! velocity
   real, intent(inout)    :: w(nx,ny,nz)                       ! velocity
   integer, intent(in)    :: it
   real, intent(in)       :: ampl
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
#endif

   integer lit,j,k,ip,l,ii
#ifdef _CUDA
   type(dim3) :: grid,tblock
#endif

   real, parameter :: inv1cs2 = 1.0/(cs2)
   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)

   integer, parameter :: icpu=3
   call cpustart()


   ip=iturb_pos
   lit=mod(it,nrturb)
   if (lit==0) lit=nrturb

   ii = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing the S_i term returned in df
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(rho, u, v, w, ip, rtmp, vel)
#endif
   do k=1,nz
   do j=1,ny
      rtmp(1,j,k)=rho(ip,j,k)
      vel(1,1,j,k)=u(ip,j,k)
      vel(2,1,j,k)=v(ip,j,k)
      vel(3,1,j,k)=w(ip,j,k)
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute equilibrium distribution for the turbine grid points
#ifdef _CUDA
      tx=1; bx=(ii+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=8; bz=(nz+tz-1)/tz
#endif

      call fequil3_block_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(dfeq1,rtmp,vel,ii,H2,H3,cxs,cys,czs,cs2,weights,inv1cs2,inv2cs4,inv6cs6,ibgk)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(j,k) SHARED(rho, ip, rtmp, vel,  uu, vv, ww, lit)
#endif
   do k=1,nz
   do j=1,ny
      rtmp(1,j,k)=rho(ip,j,k)
      vel(1,1,j,k)=vel(1,1,j,k)+ampl*uu(j,k,lit)
      vel(2,1,j,k)=vel(2,1,j,k)+ampl*vv(j,k,lit)
      vel(3,1,j,k)=vel(3,1,j,k)+ampl*ww(j,k,lit)
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute equilibrium distribution for the turbine grid points
#ifdef _CUDA
      tx=8; bx=(ii+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=8; bz=(nz+tz-1)/tz
#endif

      call fequil3_block_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(dfeq2,rtmp,vel,ii,H2,H3,cxs,cys,czs,cs2,weights,inv1cs2,inv2cs4,inv6cs6,ibgk)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(j,k,l) SHARED(turbulence_df, dfeq1,dfeq2)
#endif
   do k=1,nz
   do j=1,ny
      do l=1,nl
         turbulence_df(l,j,k)=dfeq2(l,1,j,k)-dfeq1(l,1,j,k)
      enddo
   enddo
   enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

   call cpufinish(icpu)

end subroutine
end module
