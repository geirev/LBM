module m_boundary_i_inflow_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_i_inflow_kernel(f,uvel,udir,taperj,taperk)
! Inflow/outflow boundary conditions in the i-direction.
!
!  - For positive ux, inflow is imposed at i=0 and outflow at i=nx+1.
!  - For negative ux, inflow is imposed at i=nx+1 and outflow at i=0.
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
   real, intent(in)    :: taperj(ny)
   real, intent(in)    :: taperk(nz)
   real, value         :: udir

   integer :: j,k,l
   real, parameter :: pi = acos(-1.0)
   real :: wl, cxl, cyl
   real :: uu
   real :: fghost(nl) ! Local buffer for reconstructed ghost distributions
   real :: uxdir, uydir
   real, parameter :: dir_tol = 10.0*epsilon(1.0)
   real, parameter :: switch_tol = 0.02
   real :: alphai, rholoc


!------------------ Indexing (CUDA vs CPU) -------------------------
#ifdef _CUDA
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP& PRIVATE(j,k,l,wl,cxl,cyl,uu,fghost,uxdir,uydir,alphai,rholoc) &
!$OMP& SHARED(f,uvel,udir,taperj,taperk)
   do k = 1, nz
   do j = 1, ny
#endif
      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)
      alphai = min(1.0, abs(uxdir)/switch_tol)

! Inflow at ghost plane i=0; outflow at ghost plane i=nx+1
      if (uxdir > dir_tol) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (i=0) using
         !    Krüger-type velocity condition based on interior f(1,j,k).
         !    Store temporarily in fghost(l) to avoid device aliasing issues.
         !-----------------------------------------------------------
         uu = uvel(k)*taperj(j)*taperk(k)

         rholoc = 0.0
         do l = 1, nl
            rholoc = rholoc + f(l,1,j,k)
         enddo

         ! Krüger-style correction:
         !    fghost(l) = f(l,1,j,k) - 2 w rho (c·u)/cs2
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
         ! 3) Outflow at i=nx+1: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,nx+1,j,k) = f(l,nx,j,k)
         enddo

! Inflow at ghost plane i=nx+1; outflow at ghost plane i=0
      elseif (uxdir < -dir_tol) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (i=nx+1) using
         !    Krüger-type velocity condition based on interior f(nx,j,k).
         !    Store temporarily in fghost(l) to avoid device aliasing issues.
         !-----------------------------------------------------------
         uu = uvel(k)*taperj(j)*taperk(k)

         rholoc = 0.0
         do l = 1, nl
            rholoc = rholoc + f(l,nx,j,k)
         enddo

         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))

            ! Krüger-style correction:
            !    fghost(l) = f(l,nx,j,k) - 2 w rho (c·u)/cs2
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
         ! 3) Outflow at i=0: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,0,j,k) = f(l,1,j,k)
         enddo

      else

         !-----------------------------------------------------------
         ! 3) The normal velocity is essentially zero. Use zero-gradient
         !    extrapolation at both i-boundaries.
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

