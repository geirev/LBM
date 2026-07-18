module m_boundary_j_inflow_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_j_inflow_kernel(f,uvel,udir,taperi,taperk)
! Inflow/outflow boundary conditions in the j-direction.
!
!  - For positive uy, inflow is imposed at j=0 and outflow at j=ny+1.
!  - For negative uy, inflow is imposed at j=ny+1 and outflow at j=0.
!  - The inflow uses a Krüger-type reconstruction based on the adjacent
!    interior plane and the prescribed velocity profile uvel(k).
!  - A general opposite-direction mapping is applied using bounce(l).
!  - Open outflow uses first-order zero-gradient extrapolation.

#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions
   use mod_D3Q27setup

   implicit none
   real, intent(inout) :: f     (nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: uvel  (nz)
   real, intent(in)    :: taperi(nx)
   real, intent(in)    :: taperk(nz)
   real, value         :: udir

   integer :: i,k,l
   real, parameter :: pi = acos(-1.0)
   real :: wl, cxl, cyl
   real :: uu
   real :: fghost(nl) ! Local buffer for reconstructed ghost distributions
   real :: uxdir, uydir
   real, parameter :: dir_tol = 10.0*epsilon(1.0)
   real, parameter :: switch_tol = 0.02
   real :: alphaj, rholoc

!------------------ Indexing (CUDA vs CPU) -------------------------
#ifdef _CUDA
   i = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (i < 1 .or. i > nx) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP& PRIVATE(i,k,l,wl,cxl,cyl,uu,fghost,uxdir,uydir,alphaj,rholoc) &
!$OMP& SHARED(f,uvel,udir,taperi,taperk)
   do k = 1, nz
   do i = 1, nx
#endif
      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)
      alphaj = min(1.0, abs(uydir)/switch_tol)

! Inflow at ghost plane j=0; outflow at ghost plane j=ny+1
      if (uydir > dir_tol) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (j=0) using
         !    Krüger-type velocity condition based on interior f(i,1,k).
         !    Store temporarily in fghost(l) to avoid device aliasing issues.
         !-----------------------------------------------------------
         uu = uvel(k)*taperi(i)*taperk(k)

         rholoc = 0.0
         do l = 1, nl
            rholoc = rholoc + f(l,i,1,k)
         enddo

         ! Krüger-style correction:
         !    fghost(l) = f(l,i,1,k) - 2 w rho (c·u)/cs2
         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))
            fghost(l) = f(l,i,1,k) - 2.0*wl*rholoc * (cxl*uu*uxdir + cyl*uu*uydir)/cs2
         enddo

         !-----------------------------------------------------------
         ! 2) General y-bounce mapping on ghost plane j=0
         !-----------------------------------------------------------
         do l = 1,nl
            if (cys(l) <= 0) then
               f(l,i,0,k) = alphaj*fghost(l) + (1.0-alphaj)*f(l,i,1,k)
            else
               f(l,i,0,k) = alphaj*fghost(bounce(l)) + (1.0-alphaj)*f(l,i,1,k)
            endif
         enddo

         !-----------------------------------------------------------
         ! 3) Outflow at j=ny+1: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,i,ny+1,k) = f(l,i,ny,k)
         enddo

      ! Inflow at ghost plane j=ny+1; outflow at ghost plane j=0
      elseif (uydir < -dir_tol) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (j=ny+1) using
         !    Krüger-type velocity condition based on interior f(i,ny,k).
         !    Store temporarily in fghost(l) to avoid device aliasing issues.
         !-----------------------------------------------------------
         uu = uvel(k)*taperi(i)*taperk(k)

         rholoc = 0.0
         do l = 1, nl
            rholoc = rholoc + f(l,i,ny,k)
         enddo

         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))

            ! Krüger-style correction:
            !    fghost(l) = f(l,i,ny,k) - 2 w rho (c·u)/cs2
            fghost(l) = f(l,i,ny,k) - 2.0*wl*rholoc*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
         enddo

         !-----------------------------------------------------------
         ! 2) General y-bounce mapping on ghost plane j=ny+1
         !-----------------------------------------------------------

         do l = 1,nl
            if (cys(l) >= 0) then
               f(l,i,ny+1,k) = alphaj*fghost(l) + (1.0-alphaj)*f(l,i,ny,k)
            else
               f(l,i,ny+1,k) = alphaj*fghost(bounce(l)) + (1.0-alphaj)*f(l,i,ny,k)
            endif
         enddo


         !-----------------------------------------------------------
         ! 3) Outflow at j=0: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,i,0,k) = f(l,i,1,k)
         enddo

      else

         !-----------------------------------------------------------
         ! 3) The normal velocity is essentially zero. Use zero-gradient
         !    extrapolation at both j-boundaries.
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,i,ny+1,k) = f(l,i,ny,k)
            f(l,i,0,k) = f(l,i,1,k)
         enddo

      endif

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module

