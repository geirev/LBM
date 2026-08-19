module m_fequil_difference

contains

#ifdef _CUDA
attributes(device) &
#endif
function fequil_difference(l,dens,ux,uy,uz,inv1cs2,inv6cs6,ibgk) result(dfeq)

!-----------------------------------------------------------------------
! Equilibrium difference between a lattice population l and its
! opposite population bounce(l):
!
!       dfeq = feq(l) - feq(bounce(l))
!
! Consistent with the equilibrium used in m_postcoll_kernel.
!
! For the Hermite equilibrium:
!
!   zeroth order : cancels
!   first order  : odd  -> contributes
!   second order : even -> cancels
!   third order  : odd  -> contributes if ibgk == 3
!
! Thus H2, H3, A2 and A3 do not need to be constructed explicitly.
!
! The third-order Hermite contraction is
!
!   H3 : uuu =
!
!       (c.u)^3 - 3*cs2*(c.u)*(u.u)
!
! giving
!
!   feq(l)-feq(bounce(l)) =
!
!       2*w(l)*rho * [
!
!           (c.u)*inv1cs2
!
!         + ((c.u)^3 - 3*cs2*(c.u)*(u.u))*inv6cs6
!
!       ]
!
! for ibgk == 3.
!
! For ibgk /= 3 only the first-order odd contribution remains.
!
!-----------------------------------------------------------------------

   use mod_D3Q27setup, only : weights,cxs,cys,czs,cs2

   implicit none

   integer, value :: l
   integer, value :: ibgk

   real, value :: dens
   real, value :: ux
   real, value :: uy
   real, value :: uz

   real, value :: inv1cs2
   real, value :: inv6cs6

   real :: dfeq

   real :: cu
   real :: usq
   real :: tmp


!-----------------------------------------------------------------------
! c.u
!-----------------------------------------------------------------------

   cu = real(cxs(l))*ux + &
        real(cys(l))*uy + &
        real(czs(l))*uz


!-----------------------------------------------------------------------
! First-order odd contribution.
!-----------------------------------------------------------------------

   tmp = cu*inv1cs2


!-----------------------------------------------------------------------
! Third-order odd contribution.
!
! This is exactly the contraction of your H3 tensor with
!
!       rho*u_p*u_q*u_r*inv6cs6
!
! but evaluated analytically.
!-----------------------------------------------------------------------

   if (ibgk == 3) then

      usq = ux*ux + uy*uy + uz*uz

      tmp = tmp + &
            (cu*cu*cu - 3.0*cs2*cu*usq)*inv6cs6

   endif


!-----------------------------------------------------------------------
! Difference between opposite equilibrium populations.
!-----------------------------------------------------------------------

   dfeq = 2.0*weights(l)*dens*tmp


end function fequil_difference

end module m_fequil_difference
