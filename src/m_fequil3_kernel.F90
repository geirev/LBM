! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
module m_fequil3_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine fequil3_kernel(feq, rho, u, v, w, nx, ny, nz, nl, H2, H3, cxs, cys, czs, cs2, weights,&
                             inv1cs2, inv2cs4, inv6cs6, ibgk)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value       :: nx, ny, nz, nl
   real, intent(out)    :: feq(nl,nx+2,ny+2,nz+2)
   real, intent(in)     :: rho(nx,ny,nz)
   real, intent(in)     :: u(nx,ny,nz)
   real, intent(in)     :: v(nx,ny,nz)
   real, intent(in)     :: w(nx,ny,nz)
   real, intent(in)     :: H2(3,3,nl)
   real, intent(in)     :: H3(3,3,3,nl)
   real, intent(in)     :: weights(nl)
   integer, intent(in)  :: cxs(nl)
   integer, intent(in)  :: cys(nl)
   integer, intent(in)  :: czs(nl)
   real, value          :: inv1cs2
   real, value          :: inv2cs4
   real, value          :: inv6cs6
   real, value          :: cs2
   integer, value       :: ibgk


   real :: A0_2(3,3)
   real :: A0_3(3,3,3)
   real :: vel(3)
   real :: ratio
   real :: vratio

   integer :: i, j, k, l, p, q, r, i1, j1, k1
#ifdef _CUDA
   attributes(device) :: feq
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: H2
   attributes(device) :: H3
   attributes(device) :: cxs
   attributes(device) :: cys
   attributes(device) :: czs
   attributes(device) :: weights
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
   ratio = inv6cs6 / inv2cs4
#else
   ratio = inv6cs6 / inv2cs4
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k, i1, j1, k1, vel, A0_2, A0_3, ratio, vratio) &
!$OMP          SHARED(feq, rho, u, v, w, nx, ny, nz, nl, H2, H3, cxs, cys, czs, cs2, weights, inv1cs2, inv2cs4, inv6cs6, ibgk)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      i1=i+1
      j1=j+1
      k1=k+1

! Copy u,v,w to vel(1:3)
      vel(1)=u(i,j,k)
      vel(2)=v(i,j,k)
      vel(3)=w(i,j,k)

! A0_2
      do q=1,3
      do p=1,3
         A0_2(p,q)=rho(i,j,k)*vel(p)*vel(q)*inv2cs4
      enddo
      enddo


! 1. order
      do l=1,nl
         feq(l,i1,j1,k1) = rho(i,j,k) * (cs2 + real(cxs(l))*vel(1)&
                                             + real(cys(l))*vel(2)&
                                             + real(czs(l))*vel(3)) * inv1cs2
      enddo


! 2nd order equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
      do l=1,nl
         do q=1,3
         do p=1,3
            feq(l,i1,j1,k1)=feq(l,i1,j1,k1) + H2(p,q,l)*A0_2(p,q)
         enddo
         enddo
      enddo


! Third order components
     if (ibgk == 3) then
! compute A0_3 using A0_2 to save arithmetic
         do r=1,3
         vratio = vel(r) * ratio        ! combine vel(r) with ratio to 1 value
         do q=1,3
         do p=1,3
            A0_3(p,q,r) = A0_2(p,q) * vratio
         enddo
         enddo
         enddo

!         do r=1,3
!         do q=1,3
!         do p=1,3
!            A0_3(p,q,r)=rho(i,j,k)*vel(p)*vel(q)*vel(r)*inv6cs6
!         enddo
!         enddo
!         enddo


         do l=2,nl
            do r=1,3
            do q=1,3
            do p=1,3
               feq(l,i1,j1,k1)=feq(l,i1,j1,k1) + H3(p,q,r,l)*A0_3(p,q,r)
            enddo
            enddo
            enddo
         enddo
       endif

!! !$CUF UNROLL
!! dir$ unroll
!      do l=1,nl
!         f(l,i,j,k)= weights(l)*f(l,i,j,k)
!      enddo
       feq( 1,i1,j1,k1)= weights( 1)*feq( 1,i1,j1,k1)
       feq( 2,i1,j1,k1)= weights( 2)*feq( 2,i1,j1,k1)
       feq( 3,i1,j1,k1)= weights( 3)*feq( 3,i1,j1,k1)
       feq( 4,i1,j1,k1)= weights( 4)*feq( 4,i1,j1,k1)
       feq( 5,i1,j1,k1)= weights( 5)*feq( 5,i1,j1,k1)
       feq( 6,i1,j1,k1)= weights( 6)*feq( 6,i1,j1,k1)
       feq( 7,i1,j1,k1)= weights( 7)*feq( 7,i1,j1,k1)
       feq( 8,i1,j1,k1)= weights( 8)*feq( 8,i1,j1,k1)
       feq( 9,i1,j1,k1)= weights( 9)*feq( 9,i1,j1,k1)
       feq(10,i1,j1,k1)= weights(10)*feq(10,i1,j1,k1)
       feq(11,i1,j1,k1)= weights(11)*feq(11,i1,j1,k1)
       feq(12,i1,j1,k1)= weights(12)*feq(12,i1,j1,k1)
       feq(13,i1,j1,k1)= weights(13)*feq(13,i1,j1,k1)
       feq(14,i1,j1,k1)= weights(14)*feq(14,i1,j1,k1)
       feq(15,i1,j1,k1)= weights(15)*feq(15,i1,j1,k1)
       feq(16,i1,j1,k1)= weights(16)*feq(16,i1,j1,k1)
       feq(17,i1,j1,k1)= weights(17)*feq(17,i1,j1,k1)
       feq(18,i1,j1,k1)= weights(18)*feq(18,i1,j1,k1)
       feq(19,i1,j1,k1)= weights(19)*feq(19,i1,j1,k1)
#ifndef D3Q19
       feq(20,i1,j1,k1)= weights(20)*feq(20,i1,j1,k1)
       feq(21,i1,j1,k1)= weights(21)*feq(21,i1,j1,k1)
       feq(22,i1,j1,k1)= weights(22)*feq(22,i1,j1,k1)
       feq(23,i1,j1,k1)= weights(23)*feq(23,i1,j1,k1)
       feq(24,i1,j1,k1)= weights(24)*feq(24,i1,j1,k1)
       feq(25,i1,j1,k1)= weights(25)*feq(25,i1,j1,k1)
       feq(26,i1,j1,k1)= weights(26)*feq(26,i1,j1,k1)
       feq(27,i1,j1,k1)= weights(27)*feq(27,i1,j1,k1)
#endif

#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
