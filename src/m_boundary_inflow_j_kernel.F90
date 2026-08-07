module m_boundary_inflow_j_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_inflow_j_kernel(f,uvel,udir,rho0,rho_relax,j0_is_phys,jN_is_phys)
! Inflow/outflow boundary conditions in the j-direction.
! See boundary_i_inflow_kernel for the rho_relax blending rationale.
!
!  - j0_is_phys / jN_is_phys: true only if the local j=0 / j=ny+1 ghost
!    plane is a genuine physical domain boundary on this MPI rank.
!    For a rank whose j=0 or j=ny+1 is instead an inter-rank halo
!    interface, the corresponding write is skipped entirely, leaving
!    that plane holding whatever the halo exchange last placed there
!    (real neighbor-tile data) rather than overwriting it with a
!    fabricated inflow/outflow reconstruction. Serial/single-tile runs
!    pass both as .true., reproducing the original behavior exactly.
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
   logical, value       :: j0_is_phys   ! true if local j=0 is a real domain boundary
   logical, value       :: jN_is_phys   ! true if local j=ny+1 is a real domain boundary
   real                :: uconv
   integer :: i,k,l
   real, parameter :: pi = acos(-1.0)
   real :: wl, cxl, cyl
   real :: uu
   real :: fghost(nl)
   real :: uxdir, uydir
   real, parameter :: dir_tol = 10.0*epsilon(1.0)
   real, parameter :: switch_tol = 0.02
   real :: alphaj, rholoc, rholocal
   real :: cconv, invden
!------------------ Indexing (CUDA vs CPU) -------------------------
#ifdef _CUDA
   i = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (i < 1 .or. i > nx) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP& PRIVATE(i,k,l,wl,cxl,cyl,uu,fghost,uxdir,uydir,alphaj,rholoc,rholocal,cconv,invden) &
!$OMP& SHARED(f,uvel,udir)
   do k = 1, nz
   do i = 1, nx
#endif
      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)
      uconv = uvel(k)*abs(uydir)
      alphaj = min(1.0, abs(uydir)/switch_tol)
      cconv  = min(1.0, max(0.0, uconv))
      invden = 1.0/(1.0+cconv)
! Inflow at ghost plane j=0; outflow at ghost plane j=ny+1
      if (uydir > dir_tol) then
         uu = uvel(k)
         rholocal = 0.0
         do l = 1, nl
            rholocal = rholocal + f(l,i,1,k)
         enddo
         rholoc = rho_relax*rholocal + (1.0-rho_relax)*rho0
         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))
            fghost(l) = f(l,i,1,k) - 2.0*wl*rholoc * (cxl*uu*uxdir + cyl*uu*uydir)/cs2
         enddo

         if (j0_is_phys) then
            do l = 1,nl
               if (cys(l) <= 0) then
                  f(l,i,0,k) = alphaj*fghost(l) + (1.0-alphaj)*f(l,i,1,k)
               else
                  f(l,i,0,k) = alphaj*fghost(bounce(l)) + (1.0-alphaj)*f(l,i,1,k)
               endif
            enddo
         endif

         if (jN_is_phys) then
            do l = 1, nl
               f(l,i,ny+1,k) = (f(l,i,ny+1,k) + cconv*f(l,i,ny,k)) * invden
            enddo
         endif

      ! Inflow at ghost plane j=ny+1; outflow at ghost plane j=0
      elseif (uydir < -dir_tol) then
         uu = uvel(k)
         rholocal = 0.0
         do l = 1, nl
            rholocal = rholocal + f(l,i,ny,k)
         enddo
         rholoc = rho_relax*rholocal + (1.0-rho_relax)*rho0
         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))
            fghost(l) = f(l,i,ny,k) - 2.0*wl*rholoc*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
         enddo

         if (jN_is_phys) then
            do l = 1,nl
               if (cys(l) >= 0) then
                  f(l,i,ny+1,k) = alphaj*fghost(l) + (1.0-alphaj)*f(l,i,ny,k)
               else
                  f(l,i,ny+1,k) = alphaj*fghost(bounce(l)) + (1.0-alphaj)*f(l,i,ny,k)
               endif
            enddo
         endif

         if (j0_is_phys) then
            do l = 1, nl
               f(l,i,0,k) = (f(l,i,0,k) + cconv*f(l,i,1,k)) * invden
            enddo
         endif

      else

         if (jN_is_phys) then
            do l = 1, nl
               f(l,i,ny+1,k) = f(l,i,ny,k)
            enddo
         endif
         if (j0_is_phys) then
            do l = 1, nl
               f(l,i,0,k) = f(l,i,1,k)
            enddo
         endif

      endif
#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif
end subroutine
end module
