module m_fequil3
contains

subroutine fequil3(feq, rho, u, v, w, A0_2, A0_3, vel, it)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_wtime
   use m_reg_cp_vel_kernel
   use m_fequil3_A0_kernel
   use m_fequil3_1ord_kernel
   use m_fequil3_2ord_kernel
   use m_fequil3_3ord_kernel
   use m_reg_scalef_kernel

   implicit none
   real, intent(in)      :: rho(nx,ny,nz)
   real, intent(in)      :: u(nx,ny,nz)
   real, intent(in)      :: v(nx,ny,nz)
   real, intent(in)      :: w(nx,ny,nz)
   integer, intent(in)   :: it
   real, intent(out)     :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: feq
#endif

   real, intent(out)   :: A0_2(3,3,nx,ny,nz)
   real, intent(out)   :: A0_3(3,3,3,nx,ny,nz)
   real, intent(out)   :: vel(3,nx,ny,nz)

#ifdef _CUDA
   attributes(device) :: A0_2
   attributes(device) :: A0_3
   attributes(device) :: vel
#endif

   integer :: i, j, k, l, p, q, r, ia

   real, parameter :: inv1cs2 = 1.0/(cs2)
   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
   integer, parameter :: icpu=4
   real tmp
   integer :: tx, ty, tz, bx, by, bz


   call cpustart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copy u,v,w to vel(1:3)
!@cuf istat = cudaDeviceSynchronize()
      t0 = wallclock()
#ifdef _CUDA
      tx=8; bx=(nx+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=8; bz=(nz+tz-1)/tz
#endif
      call reg_cp_vel_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(vel, u, v, w, nx, ny, nz)
!@cuf istat = cudaDeviceSynchronize()
      t1 = wallclock(); walltimelocal(1)=walltimelocal(1)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Compute A0_2, A0_3
!@cuf istat = cudaDeviceSynchronize()
      t0 = wallclock()
#ifdef _CUDA
      tx=8; bx=(nx+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=8; bz=(nz+tz-1)/tz
#endif
      call fequil3_A0_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(rho, A0_2, A0_3, vel, nx, ny, nz, inv2cs4, inv6cs6)
!@cuf istat = cudaDeviceSynchronize()
      t1 = wallclock(); walltimelocal(2)=walltimelocal(2)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  1st order equilibrium distribution

!@cuf istat = cudaDeviceSynchronize()
      t0 = wallclock()
#ifdef _CUDA
      tx=8; bx=(nx+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=8; bz=(nz+tz-1)/tz
#endif
      call fequil3_1ord_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(feq, rho, vel, nx, ny, nz, nl, cxs, cys, czs, cs2, inv1cs2)
!@cuf istat = cudaDeviceSynchronize()
      t1 = wallclock(); walltimelocal(3)=walltimelocal(3)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  2nd order equilibrium distribution

!@cuf istat = cudaDeviceSynchronize()
      t0 = wallclock()
#ifdef _CUDA
      tx=8; bx=(nx+2+tx-1)/tx
      ty=8; by=(ny+2+ty-1)/ty
      tz=8; bz=(nz+2+tz-1)/tz
#endif
      call fequil3_2ord_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(feq, H2, A0_2, nx+2, ny+2, nz+2, nl)
!@cuf istat = cudaDeviceSynchronize()
      t1 = wallclock(); walltimelocal(4)=walltimelocal(4)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  3nd order equilibrium distribution
!  the above identically recovers the BGK equilibrium, now we add third order contributions
   if (ibgk == 3) then

!@cuf istat = cudaDeviceSynchronize()
      t0 = wallclock()
#ifdef _CUDA
      tx=8; bx=(nx+tx-1)/tx
      ty=8; by=(ny+ty-1)/ty
      tz=8; bz=(nz+tz-1)/tz
#endif
      call fequil3_3ord_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(feq, H3, A0_3, nx, ny, nz, nl)
!@cuf istat = cudaDeviceSynchronize()
      t1 = wallclock(); walltimelocal(5)=walltimelocal(5)+t1-t0

   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scaling f by the weights
!@cuf istat = cudaDeviceSynchronize()
      t0 = wallclock()
#ifdef _CUDA
      tx=8; bx=(nx+2+tx-1)/tx
      ty=8; by=(ny+2+ty-1)/ty
      tz=8; bz=(nz+2+tz-1)/tz
#endif
      call reg_scalef_kernel&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(feq, weights, nx+2, ny+2, nz+2, nl)
!@cuf istat = cudaDeviceSynchronize()
   t1 = wallclock(); walltimelocal(6)=walltimelocal(6)+t1-t0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call cpufinish(icpu)
   if (it==999) then
      do j=1,6
         print '(a,i3,g13.5)','fequil3     :',j,walltimelocal(j)
      enddo
      print '(a,g13.5)',      'fequil3     :',sum(walltimelocal(1:6))
   endif


end subroutine
end module
