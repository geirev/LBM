module m_boundary_i_inflow_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_i_inflow_kernel(f,uvel,rho0,udir,jbnd,kbnd,taperj,taperk)
! Inflow / outflow boundary conditions in i-direction.
!
!  - Inflow imposed at ghost plane i=0 using a Krüger-type formula
!    based on interior values at i=1 and a prescribed velocity profile uvel(k).
!  - A general x-bounce mapping is applied that does NOT assume any
!    particular ordering of the D3Q27 directions.
!  - Outflow at i=nx+1 uses simple zero-gradient extrapolation.

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
   real, value         :: rho0
   integer, value      :: jbnd
   integer, value      :: kbnd

   integer :: i,j,k,l,m,ka
   real, parameter :: pi = 3.1415927410125732
   real :: wl, cxl, cyl
   real :: uu
   ! Local buffer for ghost distributions at i=0
   real :: fghost(nl)
   real :: uxdir, uydir
   real, parameter :: dir_tol = 1.0e-4
   real, parameter :: switch_tol = 0.02
   real :: alphai, alphaj

   i = 1   ! interior inflow plane

!------------------ Indexing (CUDA vs CPU) -------------------------
#ifdef _CUDA
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP& PRIVATE(j,k,l,m,ka,wl,cxl,cyl,uu,fghost,uxdir,uydir) &
!$OMP& SHARED(f,uvel,rho0,udir,taperj,taperk)
   do k = 1, nz
   do j = 1, ny
#endif
      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)
      alphai = min(1.0, abs(uxdir)/switch_tol)
      alphaj = min(1.0, abs(uydir)/switch_tol)

      !inflow at i=1 and outflow at i=nx
      if (uxdir > dir_tol) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (i=0) using
         !    Krüger-type velocity condition based on interior f(1,j,k).
         !    Store temporarily in fghost(l) to avoid device aliasing issues.
         !-----------------------------------------------------------
         ka = min(max(k,1), nz)
         uu = uvel(ka)*taperj(j)*taperk(k)

         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))

            ! Krüger-style correction:
            !    fghost(l) = f(l,1,j,k) - 2 w rho (c·u)/cs2
            fghost(l) = f(l,1,j,k) - 2.0 * wl * rho0 * &
                        ( cxl*uu*uxdir + &
                          cyl*uu*uydir ) / cs2
         enddo

         !-----------------------------------------------------------
         ! 2) General x-bounce mapping on ghost plane i=0
         !-----------------------------------------------------------
         do l = 1, nl
            if (cxs(l) <= 0) then
!               f(l,0,j,k) = fghost(l)
               f(l,0,j,k) = alphai*fghost(l) + (1.0-alphai)*f(l,1,j,k)
            else
               do m = 1, nl
                  if (cxs(m) == -cxs(l) .and. &
                      cys(m) == -cys(l) .and. &
                      czs(m) == -czs(l)) then
!                     f(l,0,j,k) = fghost(m)
                     f(l,0,j,k) = alphai*fghost(m) + (1.0-alphai)*f(l,1,j,k)
                     exit
                  endif
               enddo
            endif
         enddo

         !-----------------------------------------------------------
         ! 3) Outflow at i=nx+1: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,nx+1,j,k) = f(l,nx,j,k)
         enddo

      !inflow at i=nx and outflow at i=1
      elseif (uxdir < -dir_tol) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (i=nx+1) using
         !    Krüger-type velocity condition based on interior f(nx,j,k).
         !    Store temporarily in fghost(l) to avoid device aliasing issues.
         !-----------------------------------------------------------
         ka = min(max(k,1), nz)
         uu = uvel(ka)*taperj(j)*taperk(k)

         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))

            ! Krüger-style correction:
            !    fghost(l) = f(l,nx,j,k) - 2 w rho (c·u)/cs2
            fghost(l) = f(l,nx,j,k) - 2.0 * wl * rho0 * &
                        ( cxl*uu*uxdir + &
                          cyl*uu*uydir ) / cs2
         enddo

         !-----------------------------------------------------------
         ! 2) General x-bounce mapping on ghost plane i=nx+1
         !-----------------------------------------------------------
         do l = 1, nl
            if (cxs(l) >= 0) then
!               f(l,nx+1,j,k) = fghost(l)
               f(l,nx+1,j,k) = alphai*fghost(l) + (1.0-alphai)*f(l,nx,j,k)
            else
               do m = 1, nl
                  if (cxs(m) == -cxs(l) .and. &
                      cys(m) == -cys(l) .and. &
                      czs(m) == -czs(l)) then
!                     f(l,nx+1,j,k) = fghost(m)
                     f(l,nx+1,j,k) = alphai*fghost(m) + (1.0-alphai)*f(l,nx,j,k)
                     exit
                  endif
               enddo
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
         ! 3) Outflow at i=nx+1: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,nx+1,j,k) = f(l,nx,j,k)
         enddo

         !-----------------------------------------------------------
         ! 3) Outflow at i=0: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,0,j,k) = f(l,1,j,k)
         enddo
      endif

!!      !-----------------------------------------------------------
!!      ! 1) Build "raw" inflow at ghost plane (i=0) using
!!      !    Krüger-type velocity condition based on interior f(1,j,k).
!!      !    Store temporarily in fghost(l) to avoid device aliasing issues.
!!      !-----------------------------------------------------------
!!      ka = min(max(k,1), nz)
!!      uu = uvel(ka)*taperj(j)*taperk(k)
!!
!!      do l = 1, nl
!!         wl  = weights(l)
!!         cxl = real(cxs(l))
!!         cyl = real(cys(l))
!!
!!         ! Krüger-style correction:
!!         !    fghost(l) = f(l,1,j,k) - 2 w rho (c·u)/cs2
!!         fghost(l) = f(l,1,j,k) - 2.0 * wl * rho0 * &
!!                     ( cxl*uu*cos(udir*pi/180.0) + &
!!                       cyl*uu*sin(udir*pi/180.0) ) / cs2
!!      enddo
!!
!!      !-----------------------------------------------------------
!!      ! 2) General x-bounce mapping on ghost plane i=0
!!      !
!!      !    Idea (mimicking a tmp-swap logic in a robust way):
!!      !    - for directions with cxs > 0 (incoming from ghost to fluid):
!!      !         f(l,0,j,k) := fghost(l_opp)
!!      !      where l_opp has cxs = -cxs(l),  cys=-cys(l), czs=-czs(l)
!!      !      Of course we assume no vertical component to the inflow so the czs should not matter.
!!      !
!!      !    This gives "one-timestep bounce-back" behaviour for x,
!!      !    without relying on even/odd indexing or pair ordering.
!!      !-----------------------------------------------------------
!!      do l = 1, nl
!!         if (cxs(l) <= 0) then
!!            ! keep the Krüger-corrected value directly
!!            f(l,0,j,k) = fghost(l)
!!         else
!!            ! find opposite direction in x
!!            do m = 1, nl
!!               if (cxs(m) == -cxs(l) .and. &
!!                   cys(m) == -cys(l) .and. &
!!                   czs(m) == -czs(l)) then
!!                  f(l,0,j,k) = fghost(m)
!!                  exit
!!               endif
!!            enddo
!!         endif
!!      enddo
!!
!!      !-----------------------------------------------------------
!!      ! 3) Outflow at i=nx+1: zero-gradient extrapolation
!!      !-----------------------------------------------------------
!!      do l = 1, nl
!!         f(l,nx+1,j,k) = f(l,nx,j,k)
!!      enddo


#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module

