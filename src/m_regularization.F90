module m_regularization
contains

subroutine regularization(f, feq, u, v, w, A1_2, A1_3, vel, it, nt1)
   use mod_dimensions
   use mod_D3Q27setup
   use m_ablim
   use m_readinfile
   use m_wtime

#ifdef _CUDA
    use cudafor
#endif
   use m_reg_subtract_feq_kernel
   use m_reg_cp_vel_kernel
   use m_reg_A1_2_kernel
   use m_reg_A1_3_kernel
   use m_reg_scaleA1_kernel
   use m_reg_2ord_kernel
   use m_reg_3ord_kernel
   use m_reg_scalef_kernel

   implicit none
   real, intent(in)       :: u(nx,ny,nz)
   real, intent(in)       :: v(nx,ny,nz)
   real, intent(in)       :: w(nx,ny,nz)
   real, intent(in)       :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, intent(in)   :: it
   integer, intent(in)   :: nt1
#ifdef _CUDA
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: f
   attributes(device) :: feq
#endif

   real, intent(out)   :: A1_2(3,3,nx,ny,nz)
   real, intent(out)   :: A1_3(3,3,3,nx,ny,nz)
   real, intent(out)   :: vel(1:3,nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: A1_2
   attributes(device) :: A1_3
   attributes(device) :: vel
#endif

   integer :: i, j, k, l, p, q, r

   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
   integer, parameter :: icpu=5

! Define your dimensions
   integer :: tx, ty, tz, bx, by, bz

   call cpustart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing non-equilibrium distribution defined in \citet{fen21a} between Eqs (32) and (33)
#ifdef _CUDA
   tx=ntx; bx=(nx+2+tx-1)/tx
   ty=nty; by=(ny+2+ty-1)/ty
   tz=ntz; bz=(nz+2+tz-1)/tz
#endif
   call reg_subtract_feq_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, feq, nx+2, ny+2, nz+2, nl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! projecting non-equilibrium distribution on the Hermitian polynomials for regularization
   if (ihrr == 1) then

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
!! Eq (11) from  Jacob 2018 is identical to the 33a from Feng (2021)
!! Used for regularization and turbulence calculation
!              call dgemv('n', 9,27,1.0,H2, 9,f(1,i,j,k),1,0.0,A1_2,1)
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call reg_A1_2_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(A1_2, H2, f, nx, ny, nz, nl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A1_3(HRR) from \citet{fen21a}, as defined after Eq. (34). Proof by Malaspinas 2015, Appendix B
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call reg_A1_3_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(A1_2, A1_3, vel, nx, ny, nz)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! scale A1_2 and A1_3 to save flops
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call reg_scaleA1_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(A1_2, A1_3, nx, ny, nz, inv2cs4, inv6cs6)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rfneq from \citet{fen21a}, as defined in Eq. (34)
!              call dgemv('t', 9,27,inv2cs4,H2, 9,A1_2,1,0.0,f(1,i,j,k),1)
!              call dgemv('t',27,27,inv6cs6,H3,27,A1_3,1,1.0,f(1,i,j,k),1)
!              f(:,i,j,k)=weights(:)*f(:,i,j,k)
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call reg_2ord_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, A1_2, H2, nx, ny, nz, nl)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=ntz; bz=(nz+2+tz-1)/tz
#endif
      call reg_3ord_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, A1_3, H3, nx+2, ny+2, nz+2, nl)



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
        &(f, weights, nx+2, ny+2, nz+2, nl)

   endif
   call cpufinish(icpu)

end subroutine

end module
