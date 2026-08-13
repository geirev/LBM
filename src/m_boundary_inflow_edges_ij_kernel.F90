module m_boundary_inflow_edges_ij_kernel
contains
#ifdef _CUDA
attributes(global) &
#endif
subroutine boundary_inflow_edges_ij_kernel(f,uvel,udir,rho0,rho_relax)

#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions
   use mod_D3Q27setup, only : nl, weights, cxs, cys, bounce, cs2

   implicit none

   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: uvel(nz)

   real, value :: udir
   real, value :: rho0
   real, value :: rho_relax
   real        :: uconv

   real, parameter :: udir_tol = 10.0*epsilon(1.0)
   real, parameter :: pi       = acos(-1.0)
   real, parameter :: uconv_min_frac = 0.2   ! same floor fraction as face kernels

#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uvel
#endif

   integer :: k,l,m,mm
   real :: uxdir,uydir
   real :: uu,rhocorner,rholocal
   real :: momentum_correction
   real :: cconv, invden
   real :: uvel_ref   ! reference bulk speed for the floor (see note below)

#ifdef _CUDA
   l = threadIdx%x + (blockIdx%x-1)*blockDim%x
   k = threadIdx%y + (blockIdx%y-1)*blockDim%y
   if (l < 1 .or. l > nl) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                                        &
!$OMP PRIVATE(l,k,m,mm,uxdir,uydir,uu,rhocorner,rholocal,momentum_correction,cconv,invden,uvel_ref) &
!$OMP SHARED(f,uvel,udir,rho0,rho_relax)
   do k = 1,nz
   do l = 1,nl
#endif

      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)

      uu = uvel(k)

      ! Floor uconv against the bulk imposed speed so a locally small
      ! or zero uvel(k) (e.g. near a k=1/k=nz wall-adjacent profile)
      ! cannot make this corner's outflow ghost effectively frozen.
      ! maxval(uvel) is the natural bulk reference here since the
      ! corner's own uconv is already the full (unprojected) speed,
      ! unlike the face kernels where the floor guards against an
      ! oblique-angle projection going to zero.
      uvel_ref = maxval(uvel)
      uconv = max(uu, uconv_min_frac*uvel_ref)

      cconv  = min(1.0, max(0.0, uconv))
      invden = 1.0/(1.0+cconv)

      !------------------------------------------------------------
      ! Fallback / mixed-corner treatment
      !------------------------------------------------------------

      f(l,0,0,k)       = 0.5*(f(l,0,1,k)     + f(l,1,0,k))
      f(l,0,ny+1,k)    = 0.5*(f(l,0,ny,k)    + f(l,1,ny+1,k))
      f(l,nx+1,0,k)    = 0.5*(f(l,nx+1,1,k)  + f(l,nx,0,k))
      f(l,nx+1,ny+1,k) = 0.5*(f(l,nx+1,ny,k) + f(l,nx,ny+1,k))

      !------------------------------------------------------------
      ! Quadrant 1: ux>0, uy>0  -> inflow corner (0,0), ref node (1,1,k)
      !------------------------------------------------------------
      if (uxdir > udir_tol .and. uydir > udir_tol) then

         if (cxs(l) >= 0 .and. cys(l) >= 0) then
            rholocal = 0.0
            do mm = 1, nl
               rholocal = rholocal + f(mm,1,1,k)
            enddo
            rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0

            m = bounce(l)
            momentum_correction =                                     &
               2.0*weights(l)*rhocorner*uu                            &
               * (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
            f(l,0,0,k) = f(m,1,1,k) + momentum_correction
         endif

         f(l,nx+1,ny+1,k) =                                           &
            (f(l,nx+1,ny+1,k) + cconv*f(l,nx,ny,k)) * invden

      !------------------------------------------------------------
      ! Quadrant 4: ux>0, uy<0  -> inflow corner (0,ny+1), ref node (1,ny,k)
      !------------------------------------------------------------
      elseif (uxdir > udir_tol .and. uydir < -udir_tol) then

         if (cxs(l) >= 0 .and. cys(l) <= 0) then
            rholocal = 0.0
            do mm = 1, nl
               rholocal = rholocal + f(mm,1,ny,k)
            enddo
            rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0

            m = bounce(l)
            momentum_correction =                                     &
               2.0*weights(l)*rhocorner*uu                            &
               * (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
            f(l,0,ny+1,k) = f(m,1,ny,k) + momentum_correction
         endif

         f(l,nx+1,0,k) =                                              &
            (f(l,nx+1,0,k) + cconv*f(l,nx,1,k)) * invden

      !------------------------------------------------------------
      ! Quadrant 2: ux<0, uy>0  -> inflow corner (nx+1,0), ref node (nx,1,k)
      !------------------------------------------------------------
      elseif (uxdir < -udir_tol .and. uydir > udir_tol) then

         if (cxs(l) <= 0 .and. cys(l) >= 0) then
            rholocal = 0.0
            do mm = 1, nl
               rholocal = rholocal + f(mm,nx,1,k)
            enddo
            rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0

            m = bounce(l)
            momentum_correction =                                     &
               2.0*weights(l)*rhocorner*uu                            &
               * (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
            f(l,nx+1,0,k) = f(m,nx,1,k) + momentum_correction
         endif

         f(l,0,ny+1,k) =                                              &
            (f(l,0,ny+1,k) + cconv*f(l,1,ny,k)) * invden

      !------------------------------------------------------------
      ! Quadrant 3: ux<0, uy<0  -> inflow corner (nx+1,ny+1), ref node (nx,ny,k)
      !------------------------------------------------------------
      elseif (uxdir < -udir_tol .and. uydir < -udir_tol) then

         if (cxs(l) <= 0 .and. cys(l) <= 0) then
            rholocal = 0.0
            do mm = 1, nl
               rholocal = rholocal + f(mm,nx,ny,k)
            enddo
            rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0

            m = bounce(l)
            momentum_correction =                                     &
               2.0*weights(l)*rhocorner*uu                            &
               * (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
            f(l,nx+1,ny+1,k) = f(m,nx,ny,k) + momentum_correction
         endif

         f(l,0,0,k) =                                                 &
            (f(l,0,0,k) + cconv*f(l,1,1,k)) * invden

      endif

#ifndef _CUDA
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
