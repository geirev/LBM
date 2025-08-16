module m_fequil3
contains

subroutine fequil3(feq, rho, u, v, w, A0_2, A0_3, vel)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_wtime
   use m_reg_cp_vel_kernel
   use m_fequil3_A02_kernel
   use m_fequil3_A03_kernel
   use m_fequil3_1ord_kernel
   use m_fequil3_2ord_kernel
   use m_fequil3_3ord_kernel
   use m_reg_scalef_kernel

   implicit none
   real, intent(in)      :: rho(nx,ny,nz)
   real, intent(in)      :: u(nx,ny,nz)
   real, intent(in)      :: v(nx,ny,nz)
   real, intent(in)      :: w(nx,ny,nz)
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
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call reg_cp_vel_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(vel, u, v, w, nx, ny, nz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Compute A0_2
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call fequil3_A02_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(rho, A0_2, vel, nx, ny, nz, inv2cs4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  1st order equilibrium distribution

#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call fequil3_1ord_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(feq, rho, vel, nx, ny, nz, nl, cxs, cys, czs, cs2, inv1cs2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  2nd order equilibrium distribution
!!!            call dgemv('n', 27,9,inv2cs4,H2, 27,A0_2,1,1.0,feq(1,i,j,k),1)

#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call fequil3_2ord_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(feq, H2, A0_2, nx+2, ny+2, nz+2, nl)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  3nd order equilibrium distribution
!  the above identically recovers the BGK equilibrium, now we add third order contributions
   if (ibgk == 3) then
!  Compute A0_3
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call fequil3_A03_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(rho, A0_2, A0_3, vel, nx, ny, nz, inv2cs4, inv6cs6)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call fequil3_3ord_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(feq, H3, A0_3, nx+2, ny+2, nz+2, nl)

!!!!  call dgemv('n',27,27,inv6cs6,H3,27,A0_3,1,1.0,feq(1,i,j,k),1)
   endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scaling f by the weights
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=ntz; bz=(nz+2+tz-1)/tz
#endif
      call reg_scalef_kernel&
#ifdef _CUDA
          &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
          &(feq, weights, nx+2, ny+2, nz+2, nl)


   call cpufinish(icpu)

end subroutine
end module
