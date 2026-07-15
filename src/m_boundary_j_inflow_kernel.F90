module m_boundary_j_inflow_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine boundary_j_inflow_kernel(f,uvel,rho0,udir,ibnd,kbnd,taperi,taperk)
! Inflow / outflow boundary conditions in j-direction.
!
!  - Inflow imposed at ghost plane j=0 using a Krüger-type formula
!    based on interior values at j=1 and a prescribed velocity profile uvel(k).
!  - A general y-bounce mapping is applied that does NOT assume any
!    particular ordering of the D3Q27 directions.
!  - Outflow at j=ny+1 uses simple zero-gradient extrapolation.

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
   real, value         :: rho0
   integer, value      :: ibnd
   integer, value      :: kbnd

   integer :: i,j,k,l,m,ka
   real, parameter :: pi = 3.1415927410125732
   real :: wl, cxl, cyl
   real :: uu
   ! Local buffer for ghost distributions at j=0
   real :: fghost(nl)

   j = 1   ! interior inflow plane

!------------------ Indexing (CUDA vs CPU) -------------------------
#ifdef _CUDA
   i = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (i < 1 .or. i > nx) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,k,l,m,ka,wl,cxl,cyl,uu,fghost) &
!$OMP& SHARED(f,uvel,rho0,udir)
   do k = 1, nz
   do i = 1, nx
#endif
      if (sin(udir*pi/180.0) >= 0.0) then


         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (j=0) using
         !    Krüger-type velocity condition based on interior f(i,1,k).
         !    Store temporarily in fghost(l) to avoid device aliasing issues.
         !-----------------------------------------------------------
         ka = min(max(k,1), nz)
         uu = uvel(ka)*taperi(i)*taperk(k)

         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))

            ! Krüger-style correction:
            !    fghost(l) = f(l,i,1,k) - 2 w rho (c·u)/cs2
            fghost(l) = f(l,i,1,k) - 2.0 * wl * rho0 * &
                        ( cxl*uu*cos(udir*pi/180.0) + &
                          cyl*uu*sin(udir*pi/180.0) ) / cs2
         enddo

         !-----------------------------------------------------------
         ! 2) General y-bounce mapping on ghost plane j=0
         !-----------------------------------------------------------
         do l = 1, nl
            if (cys(l) <= 0) then
               f(l,i,0,k) = fghost(l)
            else
               do m = 1, nl
                  if (cxs(m) == -cxs(l) .and. &
                      cys(m) == -cys(l) .and. &
                      czs(m) == -czs(l)) then
                     f(l,i,0,k) = fghost(m)
                     exit
                  endif
               enddo
            endif
         enddo

         !-----------------------------------------------------------
         ! 3) Outflow at j=ny+1: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,i,ny+1,k) = f(l,i,ny,k)
         enddo

      elseif (sin(udir*pi/180.0) < 0.0) then

         !-----------------------------------------------------------
         ! 1) Build "raw" inflow at ghost plane (j=ny+1) using
         !    Krüger-type velocity condition based on interior f(i,ny,k).
         !    Store temporarily in fghost(l) to avoid device aliasing issues.
         !-----------------------------------------------------------
         ka = min(max(k,1), nz)
         uu = uvel(ka)*taperi(i)*taperk(k)

         do l = 1, nl
            wl  = weights(l)
            cxl = real(cxs(l))
            cyl = real(cys(l))

            ! Krüger-style correction:
            !    fghost(l) = f(l,i,ny,k) - 2 w rho (c·u)/cs2
            fghost(l) = f(l,i,ny,k) - 2.0 * wl * rho0 * &
                        ( cxl*uu*cos(udir*pi/180.0) + &
                          cyl*uu*sin(udir*pi/180.0) ) / cs2
         enddo

         !-----------------------------------------------------------
         ! 2) General y-bounce mapping on ghost plane j=ny+1
         !-----------------------------------------------------------
         do l = 1, nl
            if (cys(l) >= 0) then
               f(l,i,ny+1,k) = fghost(l)
            else
               do m = 1, nl
                  if (cxs(m) == -cxs(l) .and. &
                      cys(m) == -cys(l) .and. &
                      czs(m) == -czs(l)) then
                     f(l,i,ny+1,k) = fghost(m)
                     exit
                  endif
               enddo
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
         ! 3) Outflow at j=ny+1: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,i,ny+1,k) = f(l,i,ny,k)
         enddo

         !-----------------------------------------------------------
         ! 3) Outflow at j=0: zero-gradient extrapolation
         !-----------------------------------------------------------
         do l = 1, nl
            f(l,i,0,k) = f(l,i,1,k)
         enddo
      endif

!!      !-----------------------------------------------------------
!!      ! 1) Build "raw" inflow at ghost plane (j=0) using
!!      !    Krüger-type velocity condition based on interior f(i,1,k).
!!      !    Store temporarily in fghost(l) to avoid device aliasing issues.
!!      !-----------------------------------------------------------
!!      ka = min(max(k,1), nz)
!!      uu = uvel(ka)*taperi(i)*taperk(k)
!!
!!      do l = 1, nl
!!         wl  = weights(l)
!!         cxl = real(cxs(l))
!!         cyl = real(cys(l))
!!
!!         ! Krüger-style correction:
!!         !    fghost(l) = f(l,i,1,k) - 2 w rho (c·u)/cs2
!!         fghost(l) = f(l,i,1,k) - 2.0 * wl * rho0 * &
!!                     ( cxl*uu*cos(udir*pi/180.0) + &
!!                       cyl*uu*sin(udir*pi/180.0) ) / cs2
!!      enddo
!!
!!      !-----------------------------------------------------------
!!      ! 2) General y-bounce mapping on ghost plane j=0
!!      !
!!      !    Idea (mimicking a tmp-swap logic in a robust way):
!!      !    - for directions with cys <= 0: keep fghost(l) as is
!!      !    - for directions with cys > 0 (incoming from ghost to fluid):
!!      !         f(l,i,0,k) := fghost(l_opp)
!!      !      where l_opp has cxs = -cxs(l),  cys=-cys(l), czs=-czs(l)
!!      !      Of course we assume no vertical component to the inflow so the czs should not matter.
!!      !
!!      !    This gives "one-timestep bounce-back" behaviour for y,
!!      !    without relying on even/odd indexing or pair ordering.
!!      !-----------------------------------------------------------
!!      do l = 1, nl
!!         if (cys(l) <= 0) then
!!            ! keep the Krüger-corrected value directly
!!            f(l,i,0,k) = fghost(l)
!!         else
!!            ! find opposite direction in y
!!            do m = 1, nl
!!               if (cxs(m) == -cxs(l) .and. &
!!                   cys(m) == -cys(l) .and. &
!!                   czs(m) == -czs(l)) then
!!                  f(l,i,0,k) = fghost(m)
!!                  exit
!!               endif
!!            enddo
!!         endif
!!      enddo
!!
!!      !-----------------------------------------------------------
!!      ! 3) Outflow at j=ny+1: zero-gradient extrapolation
!!      !-----------------------------------------------------------
!!      do l = 1, nl
!!         f(l,i,ny+1,k) = f(l,i,ny,k)
!!      enddo


#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
