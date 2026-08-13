module m_boundary_inflow_j_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_inflow_j_kernel(f,uvel,udir,rho0,rho_relax,uvel_ref,j0_is_phys,jN_is_phys)
! Inflow/outflow boundary conditions in the j-direction.
! See boundary_inflow_i_kernel for the smoothstep-blend rationale.

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
   real, value         :: uvel_ref
   logical, value       :: j0_is_phys
   logical, value       :: jN_is_phys

   integer :: i,k,l
   real, parameter :: pi = acos(-1.0)
   real, parameter :: blend_width    = 0.2 !0.05
   real, parameter :: uconv_min_frac = 0.1
   real :: wl, cxl, cyl
   real :: uu
   real :: fraw0(nl), frawN(nl)
   real :: fkruger0(nl), fkrugerN(nl)
   real :: fconv0(nl), fconvN(nl)
   real :: uxdir, uydir
   real :: rholoc0, rholocN, rholocal0, rholocalN
   real :: uconv, cconv, invden
   real :: xi, t, w

!------------------ Indexing (CUDA vs CPU) -------------------------
#ifdef _CUDA
   i = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (i < 1 .or. i > nx) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP& PRIVATE(i,k,l,wl,cxl,cyl,uu,fraw0,frawN,fkruger0,fkrugerN,fconv0,fconvN, &
!$OMP&         uxdir,uydir,rholoc0,rholocN,rholocal0,rholocalN,uconv,cconv,invden,xi,t,w) &
!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,uvel_ref)
   do k = 1, nz
   do i = 1, nx
#endif
      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)
      uu    = uvel(k)

      uconv  = max(uu*abs(uydir), uconv_min_frac*uvel_ref)
      cconv  = min(1.0, max(0.0, uconv))
      invden = 1.0/(1.0+cconv)

      xi = min(1.0, max(-1.0, uydir/blend_width))
      t  = 0.5*(xi+1.0)
      w  = t*t*(3.0-2.0*t)

      rholocal0 = 0.0
      do l = 1, nl
         rholocal0 = rholocal0 + f(l,i,1,k)
      enddo
      rholoc0 = rho_relax*rholocal0 + (1.0-rho_relax)*rho0

      rholocalN = 0.0
      do l = 1, nl
         rholocalN = rholocalN + f(l,i,ny,k)
      enddo
      rholocN = rho_relax*rholocalN + (1.0-rho_relax)*rho0

      do l = 1, nl
         wl  = weights(l)
         cxl = real(cxs(l))
         cyl = real(cys(l))
         fraw0(l) = f(l,i,1,k)  - 2.0*wl*rholoc0*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
         frawN(l) = f(l,i,ny,k) - 2.0*wl*rholocN*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
      enddo

      do l = 1, nl
         if (cys(l) <= 0) then
            fkruger0(l) = fraw0(l)
         else
            fkruger0(l) = fraw0(bounce(l))
         endif

         if (cys(l) >= 0) then
            fkrugerN(l) = frawN(l)
         else
            fkrugerN(l) = frawN(bounce(l))
         endif
      enddo

      do l = 1, nl
         fconv0(l) = (f(l,i,0,k)    + cconv*f(l,i,1,k))  * invden
         fconvN(l) = (f(l,i,ny+1,k) + cconv*f(l,i,ny,k)) * invden
      enddo

      if (j0_is_phys) then
         do l = 1, nl
            f(l,i,0,k) = w*fkruger0(l) + (1.0-w)*fconv0(l)
         enddo
      endif
      if (jN_is_phys) then
         do l = 1, nl
            f(l,i,ny+1,k) = (1.0-w)*fkrugerN(l) + w*fconvN(l)
         enddo
      endif

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
