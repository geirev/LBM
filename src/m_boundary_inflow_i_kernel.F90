module m_boundary_inflow_i_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_inflow_i_kernel(f,uvel,udir,rho0,rho_relax,uvel_ref)
! Inflow/outflow boundary conditions in the i-direction.
!
!  - Both ghost planes (i=0, i=nx+1) are built every step as a smooth,
!    C1-continuous blend of a Kruger-type inflow reconstruction and a
!    convective (Orlanski-type) outflow update, weighted by a
!    smoothstep function of uxdir. This replaces a hard switch between
!    formulas at uxdir=0, which previously injected a discontinuity
!    into the ghost layer every time the flow direction rotated
!    through tangential to this face - visible downstream as a
!    perturbation seeding wake meandering.
!  - blend_width sets the angular half-width (in direction-cosine
!    units) of the transition zone straddling tangential flow. Wider
!    = smoother but blends inflow/outflow physics over a broader
!    angular range; narrower = sharper but less smoothing. Tune
!    against your rotation rate - it should be wide enough, in
!    timesteps, for pressure/momentum to equilibrate as the direction
!    sweeps through it.
!  - uconv is floored against uvel_ref (a hoisted reference speed,
!    e.g. maxval(uvel) computed once outside this kernel) so a small
!    or zero local uxdir projection never freezes the outflow ghost.

#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions
   use mod_D3Q27setup

   implicit none
   real, intent(inout) :: f     (nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: uvel  (nz)
   real, value         :: udir
   real, value         :: rho0
   real, value         :: rho_relax
   real, value         :: uvel_ref   ! hoisted bulk reference speed, e.g. maxval(uvel)

   integer :: j,k,l
   real, parameter :: pi = acos(-1.0)
   real, parameter :: blend_width     = 0.2 !0.05   ! tune; direction-cosine half-width
   real, parameter :: uconv_min_frac  = 0.1
   real :: wl, cxl, cyl
   real :: uu
   real :: fraw0(nl), frawN(nl)          ! raw Kruger correction at each plane
   real :: fkruger0(nl), fkrugerN(nl)    ! bounce-mapped Kruger ghost at each plane
   real :: fconv0(nl), fconvN(nl)        ! convective-outflow candidate at each plane
   real :: uxdir, uydir
   real :: rholoc0, rholocN, rholocal0, rholocalN
   real :: uconv, cconv, invden
   real :: xi, t, w


!------------------ Indexing (CUDA vs CPU) -------------------------
#ifdef _CUDA
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP& PRIVATE(j,k,l,wl,cxl,cyl,uu,fraw0,frawN,fkruger0,fkrugerN,fconv0,fconvN, &
!$OMP&         uxdir,uydir,rholoc0,rholocN,rholocal0,rholocalN,uconv,cconv,invden,xi,t,w) &
!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,uvel_ref)
   do k = 1, nz
   do j = 1, ny
#endif
      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)
      uu    = uvel(k)

      ! Convective speed: floored so a small/zero i-normal component
      ! never freezes whichever plane is currently outflow-weighted.
      uconv  = max(uu*abs(uxdir), uconv_min_frac*uvel_ref)
      cconv  = min(1.0, max(0.0, uconv))
      invden = 1.0/(1.0+cconv)

      ! Smoothstep blend weight: w->1 as uxdir>>0 (i=0 fully inflow,
      ! nx+1 fully outflow), w->0 as uxdir<<0 (roles swapped), and a
      ! continuous, zero-derivative-at-both-ends transition in between.
      xi = min(1.0, max(-1.0, uxdir/blend_width))
      t  = 0.5*(xi+1.0)
      w  = t*t*(3.0-2.0*t)

      !-----------------------------------------------------------
      ! Local density at each plane (relaxed blend, as before).
      !-----------------------------------------------------------
      rholocal0 = 0.0
      do l = 1, nl
         rholocal0 = rholocal0 + f(l,1,j,k)
      enddo
      rholoc0 = rho_relax*rholocal0 + (1.0-rho_relax)*rho0

      rholocalN = 0.0
      do l = 1, nl
         rholocalN = rholocalN + f(l,nx,j,k)
      enddo
      rholocN = rho_relax*rholocalN + (1.0-rho_relax)*rho0

      !-----------------------------------------------------------
      ! Raw Kruger correction at each plane (direction-general;
      ! valid for either sign of uxdir).
      !-----------------------------------------------------------
      do l = 1, nl
         wl  = weights(l)
         cxl = real(cxs(l))
         cyl = real(cys(l))
         fraw0(l) = f(l,1,j,k)  - 2.0*wl*rholoc0*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
         frawN(l) = f(l,nx,j,k) - 2.0*wl*rholocN*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
      enddo

      !-----------------------------------------------------------
      ! Bounce-mapped Kruger ghost set at each plane (same mapping
      ! logic as before, now computed unconditionally at both planes).
      !-----------------------------------------------------------
      do l = 1, nl
         if (cxs(l) <= 0) then
            fkruger0(l) = fraw0(l)
         else
            fkruger0(l) = fraw0(bounce(l))
         endif

         if (cxs(l) >= 0) then
            fkrugerN(l) = frawN(l)
         else
            fkrugerN(l) = frawN(bounce(l))
         endif
      enddo

      !-----------------------------------------------------------
      ! Convective-outflow candidate at each plane (old ghost value
      ! + freshly streamed interior neighbor, as before).
      !-----------------------------------------------------------
      do l = 1, nl
         fconv0(l) = (f(l,0,j,k)    + cconv*f(l,1,j,k))  * invden
         fconvN(l) = (f(l,nx+1,j,k) + cconv*f(l,nx,j,k)) * invden
      enddo

      !-----------------------------------------------------------
      ! Smooth blend at each plane. w=1: i=0 pure inflow, nx+1 pure
      ! outflow. w=0: roles swapped. Continuous through w=0.5.
      !-----------------------------------------------------------
      do l = 1, nl
         f(l,0,j,k)    = w*fkruger0(l)          + (1.0-w)*fconv0(l)
         f(l,nx+1,j,k) = (1.0-w)*fkrugerN(l)     + w*fconvN(l)
      enddo

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
