module m_boundary_inflow_i_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_inflow_i_kernel(f,uvel,udir,rho0,rho_relax)
! Inflow/outflow boundary conditions in the i-direction.
!
!  - Local density in the Krüger correction is a relaxed blend of the
!    true local density (summed from the interior node) and the fixed
!    reference rho0:
!       rholoc = rho_relax*rho_local + (1-rho_relax)*rho0
!    rho_relax=1 recovers pure local density (most accurate, least
!    stable); rho_relax=0 recovers the old fixed-rho0 behavior (most
!    stable, biased for oblique/diagonal inflow). Tune between.
!  - A general opposite-direction mapping is applied using bounce(l).
!  - Open outflow uses a convective (Orlanski-type) boundary condition:
!       f_b^(n+1) = (f_b^n + c*f_nb^(n+1)) / (1+c),   c = uconv*dt/dx

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
   real, value         :: rho_relax   ! blend factor in [0,1], see above
   real                :: uconv   ! convective (outflow) speed, lattice units (dt=dx=1)

   integer :: j,k,l
   real, parameter :: pi = acos(-1.0)
   real :: wl, cxl, cyl
   real :: uu
   real :: fghost(nl) ! Local buffer for reconstructed ghost distributions
   real :: uxdir, uydir
   real, parameter :: dir_tol = 10.0*epsilon(1.0)
   real, parameter :: switch_tol = 0.02
   real :: alphai, rholoc, rholocal
   real :: cconv, invden   ! convective BC coefficients


!------------------ Indexing (CUDA vs CPU) -------------------------
#ifdef _CUDA
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP& PRIVATE(j,k,l,wl,cxl,cyl,uu,fghost,uxdir,uydir,alphai,rholoc,rholocal,cconv,invden) &
!$OMP& SHARED(f,uvel,udir)
   do k = 1, nz
   do j = 1, ny
#endif
      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)
      uconv = uvel(k)*abs(uxdir)
      alphai = min(1.0, abs(uxdir)/switch_tol)

      ! Clip convective coefficient: c<0 (backflow) -> freeze (c=0),
      ! cap at 1 so the interior value never dominates outright.
      cconv  = min(1.0, max(0.0, uconv))
      invden = 1.0/(1.0+cconv)

! Inflow at ghost plane i=0; outflow at ghost plane i=nx+1
      if (uxdir > dir_tol) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (i=0) using
         !    Krüger-type velocity condition based on interior f(1,j,k).
         !    rholoc is a relaxed blend of local and reference density.
         !-----------------------------------------------------------
         uu = uvel(k)

         rholocal = 0.0
         do l = 1, nl
            rholocal = rholocal + f(l,1,j,k)
         enddo
         rholoc = rho_relax*rholocal + (1.0-rho_relax)*rho0

         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))
            fghost(l) = f(l,1,j,k) - 2.0*wl*rholoc * (cxl*uu*uxdir + cyl*uu*uydir)/cs2
         enddo

         !-----------------------------------------------------------
         ! 2) General x-bounce mapping on ghost plane i=0
         !-----------------------------------------------------------
         do l = 1,nl
            if (cxs(l) <= 0) then
               f(l,0,j,k) = alphai*fghost(l) + (1.0-alphai)*f(l,1,j,k)
            else
               f(l,0,j,k) = alphai*fghost(bounce(l)) + (1.0-alphai)*f(l,1,j,k)
            endif
         enddo

         !-----------------------------------------------------------
         ! 3) Outflow at i=nx+1: convective (Orlanski-type) BC.
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,nx+1,j,k) = (f(l,nx+1,j,k) + cconv*f(l,nx,j,k)) * invden
         enddo

! Inflow at ghost plane i=nx+1; outflow at ghost plane i=0
      elseif (uxdir < -dir_tol) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (i=nx+1) using
         !    Krüger-type velocity condition based on interior f(nx,j,k).
         !-----------------------------------------------------------
         uu = uvel(k)

         rholocal = 0.0
         do l = 1, nl
            rholocal = rholocal + f(l,nx,j,k)
         enddo
         rholoc = rho_relax*rholocal + (1.0-rho_relax)*rho0

         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))
            fghost(l) = f(l,nx,j,k) - 2.0*wl*rholoc*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
         enddo

         !-----------------------------------------------------------
         ! 2) General x-bounce mapping on ghost plane i=nx+1
         !-----------------------------------------------------------
         do l = 1,nl
            if (cxs(l) >= 0) then
               f(l,nx+1,j,k) = alphai*fghost(l) + (1.0-alphai)*f(l,nx,j,k)
            else
               f(l,nx+1,j,k) = alphai*fghost(bounce(l)) + (1.0-alphai)*f(l,nx,j,k)
            endif
         enddo

         !-----------------------------------------------------------
         ! 3) Outflow at i=0: convective (Orlanski-type) BC.
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,0,j,k) = (f(l,0,j,k) + cconv*f(l,1,j,k)) * invden
         enddo

      else

         !-----------------------------------------------------------
         ! Normal velocity ~0 along i: fall back to zero-gradient at
         ! both i-boundaries.
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,nx+1,j,k) = f(l,nx,j,k)
            f(l,0,j,k) = f(l,1,j,k)
         enddo

      endif

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
