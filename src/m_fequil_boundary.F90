module m_fequil_boundary

!=======================================================================
! Boundary equilibrium utilities.
!
! These functions are algebraically consistent with the equilibrium
! distribution used in m_postcoll_kernel.
!
! fequil_value:
!
!   Computes one lattice population of the equilibrium distribution:
!
!      feq = w*rho * [
!             1
!           + cu*inv1cs2
!           + (cu^2-cs2*u^2)*inv2cs4
!           + (cu^3-3*cs2*cu*u^2)*inv6cs6
!         ]
!
!   The third-order term is included only when ibgk == 3.
!
!
! fequil_difference:
!
!   Computes the difference between a lattice population and its
!   opposite equilibrium population:
!
!      dfeq = feq(l) - feq(bounce(l))
!
!   Because opposite lattice directions satisfy
!
!      c(bounce(l)) = -c(l),
!
!   all even-order Hermite contributions cancel:
!
!      zeroth order : cancels
!      first order  : odd  -> contributes
!      second order : even -> cancels
!      third order  : odd  -> contributes if ibgk == 3
!
!   Therefore fequil_difference can be evaluated much more efficiently
!   than by calling fequil_value twice.
!
! These functions are intended for use in the open-boundary kernels:
!
!   - fequil_difference : prescribed inflow reconstruction
!   - fequil_value      : non-equilibrium extrapolation outflow
!
!=======================================================================

contains


!=======================================================================
! Single-population equilibrium distribution
!=======================================================================

#ifdef _CUDA
attributes(device) &
#endif
function fequil_value(l,dens,ux,uy,uz, &
                      inv1cs2,inv2cs4,inv6cs6,ibgk) result(feq)

!-----------------------------------------------------------------------
! Compute one population of the equilibrium distribution.
!
! Algebraically identical to the Hermite equilibrium used in
! m_postcoll_kernel.
!
! Definitions:
!
!      cu  = c(l) . u
!      usq = u . u
!
! Zeroth + first + second order:
!
!      1
!    + cu/cs2
!    + (cu^2 - cs2*u^2)/(2*cs2^2)
!
! Third order:
!
!      [cu^3 - 3*cs2*cu*u^2]/(6*cs2^3)
!
! The inverse coefficients supplied to the function correspond to
!
!      inv1cs2 = 1/cs2
!      inv2cs4 = 1/(2*cs2^2)
!      inv6cs6 = 1/(6*cs2^3)
!
! The third-order contribution is included only for ibgk == 3.
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
   real, value :: inv2cs4
   real, value :: inv6cs6

   real :: feq

   real :: cu
   real :: usq
   real :: tmp


!-----------------------------------------------------------------------
! c.u and u.u
!-----------------------------------------------------------------------

   cu = real(cxs(l))*ux + &
        real(cys(l))*uy + &
        real(czs(l))*uz

   usq = ux*ux + uy*uy + uz*uz


!-----------------------------------------------------------------------
! Zeroth + first + second-order equilibrium contribution
!-----------------------------------------------------------------------

   tmp = 1.0                                  &
       + cu*inv1cs2                           &
       + (cu*cu-cs2*usq)*inv2cs4


!-----------------------------------------------------------------------
! Third-order equilibrium contribution
!
! H3 : uuu =
!
!       (c.u)^3 - 3*cs2*(c.u)*(u.u)
!-----------------------------------------------------------------------

   if (ibgk == 3) then

      tmp = tmp + &
            (cu*cu*cu - 3.0*cs2*cu*usq)*inv6cs6

   endif


!-----------------------------------------------------------------------
! Final equilibrium population
!-----------------------------------------------------------------------

   feq = weights(l)*dens*tmp


end function fequil_value



!=======================================================================
! Difference between opposite equilibrium populations
!=======================================================================

#ifdef _CUDA
attributes(device) &
#endif
function fequil_difference(l,dens,ux,uy,uz, &
                           inv1cs2,inv6cs6,ibgk) result(dfeq)

!-----------------------------------------------------------------------
! Equilibrium difference between lattice population l and its
! opposite population bounce(l):
!
!       dfeq = feq(l) - feq(bounce(l))
!
! This is algebraically consistent with the equilibrium used in
! m_postcoll_kernel.
!
!
! For opposite lattice directions:
!
!       c(bounce(l)) = -c(l)
!
! and therefore
!
!   zeroth order : cancels
!   first order  : odd  -> contributes
!   second order : even -> cancels
!   third order  : odd  -> contributes if ibgk == 3
!
! Consequently H2, H3, A2 and A3 do not need to be explicitly
! constructed.
!
!
! The third-order Hermite contraction is
!
!       H3 : uuu =
!
!           (c.u)^3 - 3*cs2*(c.u)*(u.u)
!
! giving
!
!       feq(l) - feq(bounce(l))
!
!         = 2*w(l)*rho * [
!
!               (c.u)*inv1cs2
!
!             + ((c.u)^3
!                - 3*cs2*(c.u)*(u.u))*inv6cs6
!
!           ]
!
! for ibgk == 3.
!
! For ibgk /= 3 only the first-order odd contribution remains.
!
!
! This optimized formulation is preferable to evaluating
!
!       fequil_value(l) - fequil_value(bounce(l))
!
! because the even terms are known analytically to cancel.
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
! First-order odd contribution
!-----------------------------------------------------------------------

   tmp = cu*inv1cs2


!-----------------------------------------------------------------------
! Third-order odd contribution
!
! This is exactly the contraction of the H3 tensor with
!
!       rho*u_p*u_q*u_r*inv6cs6
!
! evaluated analytically.
!-----------------------------------------------------------------------

   if (ibgk == 3) then

      usq = ux*ux + uy*uy + uz*uz

      tmp = tmp + &
            (cu*cu*cu - 3.0*cs2*cu*usq)*inv6cs6

   endif


!-----------------------------------------------------------------------
! Difference between opposite equilibrium populations
!-----------------------------------------------------------------------

   dfeq = 2.0*weights(l)*dens*tmp


end function fequil_difference


end module m_fequil_boundary
