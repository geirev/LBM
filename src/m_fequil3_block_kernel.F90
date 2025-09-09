! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
module m_fequil3_block_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine fequil3_block_kernel(feq, rho, velocity, nx, ny, nz, nl, H2, H3, cxs, cys, czs, cs2, weights,&
                             inv1cs2, inv2cs4, inv6cs6, ibgk)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value       :: nx, ny, nz, nl
   real, intent(out)    :: feq(nl,nx,ny,nz)
   real, intent(in)     :: rho(nx,ny,nz)
   real, intent(in)     :: velocity(3,nx,ny,nz)
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
   real :: cu
   real :: tmp

   integer :: i, j, k, l, p, q, r
#ifdef _CUDA
   attributes(device) :: feq
   attributes(device) :: rho
   attributes(device) :: velocity
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
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k, vel, cu, tmp, A0_2, A0_3,  vratio) &
!$OMP     SHARED(feq, rho, velocity, nx, ny, nz, nl, H2, H3, cxs, cys, czs, cs2, weights, inv1cs2, inv2cs4, inv6cs6, ratio, ibgk)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif

! Copy u,v,w to vel(1:3)
      vel(1)=velocity(1,i,j,k)
      vel(2)=velocity(2,i,j,k)
      vel(3)=velocity(3,i,j,k)

! 1. order
      do l=1,nl
         cu = real(cxs(l))*vel(1) + real(cys(l))*vel(2) + real(czs(l))*vel(3)
         feq(l,i,j,k) = rho(i,j,k) * weights(l) * (1.0 + cu*inv1cs2)
      enddo

! A0_2
      do q=1,3
      do p=1,3
         A0_2(p,q)=rho(i,j,k)*vel(p)*vel(q)*inv2cs4
      enddo
      enddo

! 2nd order equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
      do l=1,nl
         tmp=0.0
         do q=1,3
         do p=1,3
            tmp=tmp + H2(p,q,l)*A0_2(p,q)
         enddo
         enddo
         feq(l,i,j,k)=feq(l,i,j,k) + weights(l)*tmp
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
            tmp=0.0
            do r=1,3
            do q=1,3
            do p=1,3
               tmp=tmp + H3(p,q,r,l)*A0_3(p,q,r)
            enddo
            enddo
            enddo
            feq(l,i,j,k)=feq(l,i,j,k) + weights(l)*tmp
         enddo
       endif


#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
