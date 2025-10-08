module m_turbines_forcing
contains
subroutine turbines_forcing(rho,u,v,w,it)
! Returns the S_i stored in turbine_df and possibly updated velocities
   use mod_dimensions, only : nx,ny,nz
   use mod_nrel5mw,    only : nrchords,relm
   use mod_D3Q27setup
   use m_readinfile,   only : turbrpm,p2l,ipos,jpos,kpos,nturbines,ibgk
   use m_actuatorline
   !use m_actuatorline_kernel
   use m_turbines_init
   use m_turbines_compute_force
   use m_turbines_compute_force_kernel
   use m_turbines_print_force
   use m_fequil3_block_kernel
   use m_turbines_forcing_kernel_A
   use m_turbines_forcing_kernel_B
   use m_turbines_forcing_kernel_C
#ifdef _CUDA
   use cudafor
   use iso_c_binding, only: c_sizeof
#endif
   use m_wtime
   implicit none
   integer, intent(in) :: it
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
   integer ip, jp, kp, n

   real, parameter :: inv1cs2 = 1.0/(cs2)
   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)

!   real :: force_h(0:ieps,ny,nz,3)        ! work array for computing the turbine force
   real :: forceT_h(nrchords,3)             ! Tangential (driving force) along the blade
   real :: forceN_h(nrchords,3)             ! Nornal (drag force) along the blade

   real :: forceT(nrchords,3)             ! Tangential (driving force) along the blade
   real :: forceN(nrchords,3)             ! Nornal (drag force) along the blade
   real :: relm_d(nrchords)                 ! relm on device

#ifdef _CUDA
   attributes(device) :: forceT
   attributes(device) :: forceN
   attributes(device) :: relm_d
#endif



   real :: u_h(ny, nz)
   real :: v_h(ny, nz)
   real :: w_h(ny, nz)
   real :: r_h(ny, nz)

   real, save :: dtheta=0.0
   real rps
   integer j,k
   integer, parameter :: ii = 2*ieps+1

#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif


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


   do n=1,nturbines
      ip=ipos(n)
      jp=jpos(n)
      kp=kpos(n)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! My implementation of the actuator line method by Sørensen 2002 computing the force from all the turbines
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
      call actuatorline(forceN_h,forceT_h,ny,nz,jp,kp,theta,iradius,u_h,v_h,w_h,r_h)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Establising the forces from the blades on the fluid
      forceN=forceN_h
      forceT=forceT_h
      relm_d=relm
      force=0.0
#ifdef _CUDA
      tx=32; bx=(ieps+1+tx-1)/tx
      ty=8;  by=(ny+ty-1)/ty
      tz=1;  bz=(nz+tz-1)/tz
#endif
      call turbines_compute_force_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(force,forceN,forceT,theta,jp,kp,iradius,relm_d)

#ifdef _CUDA
      istat = cudaDeviceSynchronize()
#endif

!!      call turbines_compute_force(force_h,forceN_h,forceT_h,theta,jp,kp,iradius)
!!      force(0:ieps,1:ny,1:nz,1:3)=force_h(0:ieps,1:ny,1:nz,1:3)
!      if (it<2) then
!         force_h=force
!         call turbines_print_force(force_h,it)
!      endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kupershtokh forcing
! Computing the force induced velocity increments in the circular 3D rotor plane (F/rho)
#ifdef _CUDA
      tx=32; bx=(ii+tx-1)/tx
      ty=8;  by=(ny+ty-1)/ty
      tz=1;  bz=(nz+tz-1)/tz
#endif
      call turbines_forcing_kernel_A&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(u,v,w,rho,rtmp,vel,du,dv,dw,force,iradius,ip)
#ifdef _CUDA
      istat = cudaDeviceSynchronize()
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute equilibrium distribution for the turbine grid points
#ifdef _CUDA
      tx=32; bx=(ii+tx-1)/tx
      ty=8;  by=(ny+ty-1)/ty
      tz=1;  bz=(nz+tz-1)/tz
#endif
      call fequil3_block_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(dfeq1,rtmp,vel,ii,H2,H3,cxs,cys,czs,cs2,weights,inv1cs2,inv2cs4,inv6cs6,ibgk)
#ifdef _CUDA
      istat = cudaDeviceSynchronize()
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copy macro velocities plus actuator line forcing and density to vel and rtmp
#ifdef _CUDA
      tx=32; bx=(ii+tx-1)/tx
      ty=8;  by=(ny+ty-1)/ty
      tz=1;  bz=(nz+tz-1)/tz
#endif
      call turbines_forcing_kernel_B&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(u,v,w,rho,rtmp,vel,du,dv,dw,ip)
#ifdef _CUDA
      istat = cudaDeviceSynchronize()
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute equilibrium distribution for the turbine grid points with du forcing added
#ifdef _CUDA
      tx=32; bx=(ii+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=1; bz=(nz+tz-1)/tz
#endif

      call fequil3_block_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(dfeq2,rtmp,vel,ii,H2,H3,cxs,cys,czs,cs2,weights,inv1cs2,inv2cs4,inv6cs6,ibgk)
#ifdef _CUDA
      istat = cudaDeviceSynchronize()
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute final turbine_df to be used in applyturbines
#ifdef _CUDA
      tx=512; bx=(ii*nl*ny*nz+tx-1)/tx
      ty=1;   by=1
      tz=1;   bz=1
#endif
      call turbines_forcing_kernel_C&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(turbine_df,dfeq1,dfeq2,nturbines,n)
#ifdef _CUDA
      istat = cudaDeviceSynchronize()
#endif

   enddo
   call cpufinish(icpu)

end subroutine
end module
