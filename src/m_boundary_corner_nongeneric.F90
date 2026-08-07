module m_boundary_corner_nongeneric
contains
subroutine boundary_corner_nongeneric(f,uvel,udir,rho0,rho_relax,ibnd,jbnd,kbnd, &
                                       j0_is_phys,jN_is_phys)
! Three-way vertex treatment for the case where kbnd is closed and at
! least one of ibnd/jbnd is open (==1).
!
!  - j0_is_phys / jN_is_phys: true only if the local j=0 / j=ny+1 plane
!    is a genuine physical domain boundary on this MPI rank. Gates every
!    vertex write in Case A and Case C that lies on that plane, so an
!    inter-rank interface vertex is left for the halo exchange to fill
!    rather than being overwritten with a fabricated open-boundary
!    reconstruction. Case B (i open, j closed) is unaffected: j is a
!    closed wall there, not an MPI-decomposed open plane.

   use mod_dimensions
   use mod_D3Q27setup, only : nl, weights, cxs, cys, bounce, cs2
   implicit none
   real,    intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in)    :: uvel(nz)
   real,    value         :: udir
   real,    value         :: rho0
   real,    value         :: rho_relax
   integer, intent(in)    :: ibnd,jbnd,kbnd
   logical, intent(in)    :: j0_is_phys, jN_is_phys

#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uvel
#endif

   real, parameter :: pi       = acos(-1.0)
   real, parameter :: udir_tol = 10.0*epsilon(1.0)

   integer :: l,m,mm,kc,kref,jc,jref,ic,iref
   real    :: uxdir,uydir,uu,rholocal,rhocorner,momentum_correction

   ! Periodic intersections are finalized later by periodic sweeps.
   if (ibnd==0 .or. jbnd==0 .or. kbnd==0) return

   ! This routine only covers kbnd closed; kbnd==1 is unsupported
   ! elsewhere in the driver (error-stopped at entry to boundarycond).
   if (kbnd <= 10) return

   uxdir = cos(udir*pi/180.0)
   uydir = sin(udir*pi/180.0)

   !==================================================================
   ! Case A: i and j both open. Four vertices, two k-planes.
   !==================================================================
   if (ibnd==1 .and. jbnd==1) then

      do mm = 1,2
         if (mm==1) then
            kc = 0;    kref = 1
         else
            kc = nz+1; kref = nz
         endif

         uu = uvel(kref)

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
         do l = 1,nl

            ! Fallback average for populations not covered by the
            ! true-inflow override below (mirrors the ij-edge kernel).
            ! Each corner's write is gated by whether its j-plane is a
            ! genuine physical boundary on this rank.
            if (j0_is_phys) then
               f(l,0,0,kc)    = 0.5*(f(l,0,1,kc)    + f(l,1,0,kc))
               f(l,nx+1,0,kc) = 0.5*(f(l,nx+1,1,kc) + f(l,nx,0,kc))
            endif
            if (jN_is_phys) then
               f(l,0,ny+1,kc)    = 0.5*(f(l,0,ny,kc)    + f(l,1,ny+1,kc))
               f(l,nx+1,ny+1,kc) = 0.5*(f(l,nx+1,ny,kc) + f(l,nx,ny+1,kc))
            endif

            if (uxdir > udir_tol .and. uydir > udir_tol) then
               if (j0_is_phys .and. cxs(l) >= 0 .and. cys(l) >= 0) then
                  rholocal = sum(f(:,1,1,kref))
                  rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
                  m = bounce(l)
                  momentum_correction = 2.0*weights(l)*rhocorner*uu * &
                       (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
                  f(l,0,0,kc) = f(m,1,1,kref) + momentum_correction
               endif

            elseif (uxdir > udir_tol .and. uydir < -udir_tol) then
               if (jN_is_phys .and. cxs(l) >= 0 .and. cys(l) <= 0) then
                  rholocal = sum(f(:,1,ny,kref))
                  rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
                  m = bounce(l)
                  momentum_correction = 2.0*weights(l)*rhocorner*uu * &
                       (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
                  f(l,0,ny+1,kc) = f(m,1,ny,kref) + momentum_correction
               endif

            elseif (uxdir < -udir_tol .and. uydir > udir_tol) then
               if (j0_is_phys .and. cxs(l) <= 0 .and. cys(l) >= 0) then
                  rholocal = sum(f(:,nx,1,kref))
                  rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
                  m = bounce(l)
                  momentum_correction = 2.0*weights(l)*rhocorner*uu * &
                       (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
                  f(l,nx+1,0,kc) = f(m,nx,1,kref) + momentum_correction
               endif

            elseif (uxdir < -udir_tol .and. uydir < -udir_tol) then
               if (jN_is_phys .and. cxs(l) <= 0 .and. cys(l) <= 0) then
                  rholocal = sum(f(:,nx,ny,kref))
                  rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
                  m = bounce(l)
                  momentum_correction = 2.0*weights(l)*rhocorner*uu * &
                       (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
                  f(l,nx+1,ny+1,kc) = f(m,nx,ny,kref) + momentum_correction
               endif
            endif

         enddo
      enddo

   !==================================================================
   ! Case B: i open, j closed. Four vertices (ic x kc), reference the
   ! nearest real interior fluid node (iref, jref, kref).
   ! Not j0/jN-gated: jbnd>10 here, j is a closed wall on the whole
   ! domain, not an MPI-decomposed open plane.
   !==================================================================
   elseif (ibnd==1 .and. jbnd>10) then

      do mm = 1,4
         select case(mm)
         case(1); ic=0;    iref=1;  jc=0;    jref=1
         case(2); ic=0;    iref=1;  jc=ny+1; jref=ny
         case(3); ic=nx+1; iref=nx; jc=0;    jref=1
         case(4); ic=nx+1; iref=nx; jc=ny+1; jref=ny
         end select

         do kc = 0, nz+1, nz+1
            kref = merge(1,nz,kc==0)
            uu   = uvel(kref)

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
            do l = 1,nl
               if ( (ic==0    .and. cxs(l) >= 0) .or. &
                    (ic==nx+1 .and. cxs(l) <= 0) ) then
                  rholocal = sum(f(:,iref,jref,kref))
                  rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
                  m = bounce(l)
                  momentum_correction = 2.0*weights(l)*rhocorner*uu * &
                       (real(cxs(l))*uxdir)/cs2
                  f(l,ic,jc,kc) = f(m,iref,jref,kref) + momentum_correction
               else
                  f(l,ic,jc,kc) = f(l,iref,jref,kref)
               endif
            enddo
         enddo
      enddo

   !==================================================================
   ! Case C: j open, i closed. Mirror of Case B.
   !==================================================================
   elseif (jbnd==1 .and. ibnd>10) then

      do mm = 1,4
         select case(mm)
         case(1); jc=0;    jref=1;  ic=0;    iref=1
         case(2); jc=0;    jref=1;  ic=nx+1; iref=nx
         case(3); jc=ny+1; jref=ny; ic=0;    iref=1
         case(4); jc=ny+1; jref=ny; ic=nx+1; iref=nx
         end select

         if ( (jc==0    .and. .not. j0_is_phys) .or. &
              (jc==ny+1 .and. .not. jN_is_phys) ) cycle

         do kc = 0, nz+1, nz+1
            kref = merge(1,nz,kc==0)
            uu   = uvel(kref)

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
            do l = 1,nl
               if ( (jc==0    .and. cys(l) >= 0) .or. &
                    (jc==ny+1 .and. cys(l) <= 0) ) then
                  rholocal = sum(f(:,iref,jref,kref))
                  rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
                  m = bounce(l)
                  momentum_correction = 2.0*weights(l)*rhocorner*uu * &
                       (real(cys(l))*uydir)/cs2
                  f(l,ic,jc,kc) = f(m,iref,jref,kref) + momentum_correction
               else
                  f(l,ic,jc,kc) = f(l,iref,jref,kref)
               endif
            enddo
         enddo
      enddo

   endif

end subroutine
end module
