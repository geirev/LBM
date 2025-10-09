! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
module m_fequil3_block_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine fequil3_block_kernel(feq, rho, velocity, ii, H2, H3, cxs, cys, czs, cs2, weights, inv1cs2, inv2cs4, inv6cs6, ibgk)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : ny,nz
   use mod_D3Q27setup, only : nl
   implicit none
   integer, value       :: ii
   real, intent(out), contiguous    :: feq(nl,ii,ny,nz)
   real, intent(in), contiguous     :: rho(ii,ny,nz)
   real, intent(in), contiguous     :: velocity(3,ii,ny,nz)
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
   real :: tmp0,tmp1,tmp2,tmp3,rr


   integer :: i, j, k, l, p, q, r
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > ii .or. j > ny .or. k > nz) return
   ratio = inv6cs6 / inv2cs4
#else
   ratio = inv6cs6 / inv2cs4
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k, vel, cu, tmp, A0_2, A0_3,  vratio) &
!$OMP     SHARED(feq, rho, velocity, ii, H2, H3, cxs, cys, czs, cs2, weights, inv1cs2, inv2cs4, inv6cs6, ratio, ibgk)
   do k=1,nz
   do j=1,ny
   do i=1,ii
#endif

      vel(1)=velocity(1,i,j,k)
      vel(2)=velocity(2,i,j,k)
      vel(3)=velocity(3,i,j,k)
      rr=rho(i,j,k)

      do q=1,3
      do p=1,3
         A0_2(p,q)=rho(i,j,k)*vel(p)*vel(q)*inv2cs4
      enddo
      enddo

     if (ibgk == 3) then
         do r=1,3
         vratio = vel(r) * ratio
         do q=1,3
         do p=1,3
            A0_3(p,q,r) = A0_2(p,q) * vratio
         enddo
         enddo
         enddo
      endif

      do l=1,nl
         cu = real(cxs(l))*vel(1) + real(cys(l))*vel(2) + real(czs(l))*vel(3)
         tmp1 =  rho(i,j,k) * (1.0 + cu*inv1cs2)

         tmp2=0.0
         do q=1,3
         do p=1,3
            tmp2=tmp2 + H2(p,q,l)*A0_2(p,q)
         enddo
         enddo

         tmp3=0.0
         if (ibgk == 3) then
            do r=1,3
            do q=1,3
            do p=1,3
               tmp3=tmp3 + H3(p,q,r,l)*A0_3(p,q,r)
            enddo
            enddo
            enddo
         endif
         feq(l,i,j,k)= weights(l)*(tmp1 + tmp2 + tmp3)

      enddo

#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
