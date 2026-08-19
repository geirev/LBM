module m_fequil_boundary

!=======================================================================
! Single-population equilibrium distribution.
!
! Algebraically identical to the equilibrium used in postcoll_kernel:
!
!   feq = w*rho * [
!          1
!        + cu*inv1cs2
!        + (cu^2-cs2*u^2)*inv2cs4
!        + (cu^3-3*cs2*cu*u^2)*inv6cs6
!      ]
!
! The third-order term is included only for ibgk == 3.
!=======================================================================

contains


#ifdef _CUDA
attributes(device) &
#endif
function fequil_value(l,dens,ux,uy,uz, &
                      inv1cs2,inv2cs4,inv6cs6,ibgk) result(feq)

   use mod_D3Q27setup, only : weights,cxs,cys,czs,cs2

   implicit none

   integer, value :: l
   integer, value :: ibgk

   real, value :: dens
   real, value :: ux,uy,uz

   real, value :: inv1cs2
   real, value :: inv2cs4
   real, value :: inv6cs6

   real :: feq
   real :: cu
   real :: usq
   real :: tmp


   cu = real(cxs(l))*ux + &
        real(cys(l))*uy + &
        real(czs(l))*uz

   usq = ux*ux + uy*uy + uz*uz


   ! Zeroth + first + second-order equilibrium
   tmp = 1.0                                  &
       + cu*inv1cs2                           &
       + (cu*cu-cs2*usq)*inv2cs4


   ! Third-order equilibrium
   if (ibgk == 3) then

      tmp = tmp + &
            (cu*cu*cu - 3.0*cs2*cu*usq)*inv6cs6

   endif


   feq = weights(l)*dens*tmp


end function fequil_value



#ifdef _CUDA
attributes(host,device) &
#endif
function fequil_difference(l,dens,ux,uy,uz, &
                           inv1cs2,inv6cs6,ibgk) result(dfeq)

!-----------------------------------------------------------------------
! Difference between opposite equilibrium populations:
!
!       dfeq = feq(l) - feq(bounce(l))
!
! Only the odd Hermite terms survive.
!-----------------------------------------------------------------------

   use mod_D3Q27setup, only : weights,cxs,cys,czs,cs2

   implicit none

   integer, value :: l
   integer, value :: ibgk

   real, value :: dens
   real, value :: ux,uy,uz

   real, value :: inv1cs2
   real, value :: inv6cs6

   real :: dfeq
   real :: cu
   real :: usq
   real :: tmp


   cu = real(cxs(l))*ux + &
        real(cys(l))*uy + &
        real(czs(l))*uz


   tmp = cu*inv1cs2


   if (ibgk == 3) then

      usq = ux*ux + uy*uy + uz*uz

      tmp = tmp + &
            (cu*cu*cu - 3.0*cs2*cu*usq)*inv6cs6

   endif


   dfeq = 2.0*weights(l)*dens*tmp


end function fequil_difference


end module m_fequil_boundary
