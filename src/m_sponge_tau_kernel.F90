module m_sponge_tau_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine sponge_tau_kernel(tau,tau_max,nsponge_i,nsponge_j, &
                              ibnd,jbnd,j0_is_phys,jN_is_phys)
! Ramps the existing (turbulence-scheme-computed) relaxation-time
! field UP toward tau_max as nodes approach an OPEN boundary, using a
! smoothstep profile (zero derivative at both ends of the ramp).
!
! tau is intent(inout): each node's incoming value (from your
! turbulence model) is treated as the local floor and blended toward
! tau_max by the sponge weight s:
!    tau_new = tau_old + (tau_max - tau_old)*s
! s=0 in the bulk leaves tau untouched; s=1 at the open boundary pins
! tau to tau_max regardless of what the turbulence scheme computed
! there. This only increases tau (assuming tau_max exceeds whatever
! the turbulence model produces near the boundary, which should
! normally hold) - it never lowers a turbulence-model value.
!
! Only open boundaries get a sponge, with the same MPI rank-gating as
! before:
!   - i=1 and i=nx sponged whenever ibnd==1 (i never MPI-decomposed).
!   - j=1 sponged only when jbnd==1 .and. j0_is_phys (south-most rank
!     only); j=ny only when jbnd==1 .and. jN_is_phys (north-most rank
!     only). Interior ranks get no j-sponge.
!   - Closed boundaries are not sponged here.
!
! Where an i-sponge and j-sponge overlap (corner-adjacent region),
! tau uses the MAX of the two ramp weights, so the corner - the
! consistently worst-behaved region in this investigation - gets at
! least as much damping as either face alone.
!
! Call this EVERY timestep, after your turbulence scheme updates tau
! and before the collision step - unlike the earlier fixed-tau_phys
! version, this is no longer a one-time initialization call, since it
! must see and modify the current step's turbulence-computed field.

#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx, ny, nz

   implicit none
   real,    intent(inout) :: tau(0:nx+1,0:ny+1,0:nz+1)
   real,    value          :: tau_max       ! ceiling reached at the open boundary
   integer, value          :: nsponge_i     ! sponge depth in grid points, i-direction
   integer, value          :: nsponge_j     ! sponge depth in grid points, j-direction
   integer, value          :: ibnd, jbnd
   logical, value           :: j0_is_phys
   logical, value           :: jN_is_phys

#ifdef _CUDA
   attributes(device) :: tau
#endif

   integer :: i,j,k
   real :: xi_i_lo, xi_i_hi, xi_j_lo, xi_j_hi
   real :: s_i, s_j, s
   real :: tau_old

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x-1)*blockDim%x
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (i < 1 .or. i > nx) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(NONE) &
!$OMP& PRIVATE(i,j,k,xi_i_lo,xi_i_hi,xi_j_lo,xi_j_hi,s_i,s_j,s,tau_old) &
!$OMP& SHARED(tau,tau_max,nsponge_i,nsponge_j,ibnd,jbnd,j0_is_phys,jN_is_phys)
   do k = 1, nz
   do j = 1, ny
   do i = 1, nx
#endif

      !-----------------------------------------------------------
      ! i-direction sponge weight (0 in the bulk, ->1 at i=1 or i=nx)
      !-----------------------------------------------------------
      s_i = 0.0
      if (ibnd == 1) then
         xi_i_lo = 0.0
         xi_i_hi = 0.0
         if (nsponge_i > 0) then
            xi_i_lo = max(0.0, min(1.0, real(nsponge_i - (i-1))  / real(nsponge_i)))
            xi_i_hi = max(0.0, min(1.0, real(nsponge_i - (nx-i)) / real(nsponge_i)))
         endif
         s_i = max(smoothstep(xi_i_lo), smoothstep(xi_i_hi))
      endif

      !-----------------------------------------------------------
      ! j-direction sponge weight, gated per-rank by j0_is_phys/jN_is_phys
      !-----------------------------------------------------------
      s_j = 0.0
      if (jbnd == 1) then
         xi_j_lo = 0.0
         xi_j_hi = 0.0
         if (nsponge_j > 0) then
            if (j0_is_phys) xi_j_lo = max(0.0, min(1.0, real(nsponge_j - (j-1))  / real(nsponge_j)))
            if (jN_is_phys) xi_j_hi = max(0.0, min(1.0, real(nsponge_j - (ny-j)) / real(nsponge_j)))
         endif
         s_j = max(smoothstep(xi_j_lo), smoothstep(xi_j_hi))
      endif

      !-----------------------------------------------------------
      ! Combine: corner region gets the stronger of the two ramps.
      !-----------------------------------------------------------
      s = max(s_i, s_j)

      !-----------------------------------------------------------
      ! Ramp the EXISTING local tau up toward tau_max; s=0 leaves it
      ! untouched, s=1 pins it to tau_max.
      !-----------------------------------------------------------
      if (s > 0.0) then
         tau_old      = tau(i,j,k)
         tau(i,j,k)   = tau_old + (tau_max - tau_old)*s
      endif

#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
end subroutine

#ifdef _CUDA
   attributes(device) &
#endif
   pure function smoothstep(xi) result(sres)
      real, intent(in) :: xi
      real :: sres
      sres = xi*xi*(3.0 - 2.0*xi)
   end function smoothstep
end module
