module m_boundary_inflow_j_kernel
contains




!=======================================================================
! j-direction open boundary kernel
!=======================================================================

#ifdef _CUDA
attributes(global) &
#endif
subroutine boundary_inflow_j_kernel(f,uvel,udir,rho0,rho_relax, &
                                    inv1cs2,inv2cs4,inv6cs6,ibgk, &
                                    j0_is_phys,jN_is_phys)

!-----------------------------------------------------------------------
! General inflow / outflow boundary condition in j direction.
!
!
! INFLOW
! ======
!
! Incoming populations use
!
!       f_l = f_bounce(l) + feq_l - feq_bounce(l)
!
! consistently with postcoll_kernel.
!
!
! OUTFLOW
! =======
!
! Incoming populations from an outlet ghost plane use
!
!     fghost =
!
!         feq(rho0,u_int)
!
!       + [f_int - feq(rho_int,u_int)]
!
! Thus the equilibrium density is anchored to rho0 while the local
! velocity and non-equilibrium wake content are extrapolated outward.
!
!
! LOWER j FACE:
!
!   cy > 0 populations enter from j=0.
!
!   uy > 0 : inflow
!   uy < 0 : outflow
!
!
! UPPER j FACE:
!
!   cy < 0 populations enter from j=ny+1.
!
!   uy < 0 : inflow
!   uy > 0 : outflow
!
!
! Around uy=0, both reconstructions smoothly approach zero-gradient.
!
! j0_is_phys / jN_is_phys protect internal MPI j interfaces.
!-----------------------------------------------------------------------

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions

   use mod_D3Q27setup, only : nl,cxs,cys,czs,bounce

   use m_fequil_difference, only : fequil_difference
   use m_fequil_boundary, only : fequil_value

   implicit none


!=======================================================================
! Arguments
!=======================================================================

   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)

   real, intent(in) :: uvel(nz)

   real, value :: udir
   real, value :: rho0
   real, value :: rho_relax

   real, value :: inv1cs2
   real, value :: inv2cs4
   real, value :: inv6cs6

   integer, value :: ibgk

   logical, value :: j0_is_phys
   logical, value :: jN_is_phys


!=======================================================================
! Local variables
!=======================================================================

   integer :: i,k,l,m

   real, parameter :: pi = acos(-1.0)

   real, parameter :: blend_width = 0.05


!-----------------------------------------------------------------------
! Prescribed inflow velocity
!-----------------------------------------------------------------------

   real :: uxdir,uydir
   real :: uu
   real :: ux,uy


!-----------------------------------------------------------------------
! Local density and momentum at adjacent interior nodes
!-----------------------------------------------------------------------

   real :: rholocal0
   real :: rholocalN

   real :: rholoc0
   real :: rholocN

   real :: momx0,momy0,momz0
   real :: momxN,momyN,momzN

   real :: uloc0,vloc0,wloc0
   real :: ulocN,vlocN,wlocN


!-----------------------------------------------------------------------
! Smooth weights
!-----------------------------------------------------------------------

   real :: s

   real :: win0
   real :: winN

   real :: wout0
   real :: woutN


!-----------------------------------------------------------------------
! Population values
!-----------------------------------------------------------------------

   real :: fint0,fintN

   real :: fin0,finN
   real :: fout0,foutN

   real :: feq_int
   real :: feq_out

   real :: dfeq


!=======================================================================
! CUDA / CPU indexing
!=======================================================================

#ifdef _CUDA

   i = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z

   if (i < 1 .or. i > nx) return
   if (k < 1 .or. k > nz) return

#else

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                         &
!$OMP& PRIVATE(i,k,l,m,uxdir,uydir,uu,ux,uy,                       &
!$OMP&         rholocal0,rholocalN,rholoc0,rholocN,                &
!$OMP&         momx0,momy0,momz0,momxN,momyN,momzN,                &
!$OMP&         uloc0,vloc0,wloc0,ulocN,vlocN,wlocN,                &
!$OMP&         s,win0,winN,wout0,woutN,                            &
!$OMP&         fint0,fintN,fin0,finN,fout0,foutN,                  &
!$OMP&         feq_int,feq_out,dfeq)                               &
!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,                          &
!$OMP&        inv1cs2,inv2cs4,inv6cs6,ibgk,                       &
!$OMP&        j0_is_phys,jN_is_phys)

   do k = 1,nz
   do i = 1,nx

#endif


!=======================================================================
! Prescribed undisturbed inflow velocity
!=======================================================================

      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)

      uu = uvel(k)

      ux = uu*uxdir
      uy = uu*uydir


!=======================================================================
! Smooth inflow/outflow activation weights
!
! LOWER j=0:
!
!     win0  -> 1 for uy > 0
!     wout0 -> 1 for uy < 0
!
! UPPER j=ny+1:
!
!     winN  -> 1 for uy < 0
!     woutN -> 1 for uy > 0
!=======================================================================

      s = min(1.0,max(0.0, uydir/blend_width))
      win0  = s*s*(3.0-2.0*s)
      woutN = win0

      s = min(1.0,max(0.0,-uydir/blend_width))
      winN  = s*s*(3.0-2.0*s)
      wout0 = winN


!=======================================================================
! Local density and momentum at adjacent interior nodes
!=======================================================================

      rholocal0 = 0.0
      momx0     = 0.0
      momy0     = 0.0
      momz0     = 0.0

      rholocalN = 0.0
      momxN     = 0.0
      momyN     = 0.0
      momzN     = 0.0


      do l = 1,nl

         rholocal0 = rholocal0 + f(l,i,1,k)

         momx0 = momx0 + real(cxs(l))*f(l,i,1,k)
         momy0 = momy0 + real(cys(l))*f(l,i,1,k)
         momz0 = momz0 + real(czs(l))*f(l,i,1,k)


         rholocalN = rholocalN + f(l,i,ny,k)

         momxN = momxN + real(cxs(l))*f(l,i,ny,k)
         momyN = momyN + real(cys(l))*f(l,i,ny,k)
         momzN = momzN + real(czs(l))*f(l,i,ny,k)

      enddo


!=======================================================================
! Density used by prescribed inflow reconstruction
!=======================================================================

      rholoc0 = rho_relax*rholocal0 + &
                (1.0-rho_relax)*rho0

      rholocN = rho_relax*rholocalN + &
                (1.0-rho_relax)*rho0


!=======================================================================
! Local velocity used by NEE outflow reconstruction
!=======================================================================

      if (rholocal0 > tiny(1.0)) then

         uloc0 = momx0/rholocal0
         vloc0 = momy0/rholocal0
         wloc0 = momz0/rholocal0

      else

         uloc0 = ux
         vloc0 = uy
         wloc0 = 0.0

      endif


      if (rholocalN > tiny(1.0)) then

         ulocN = momxN/rholocalN
         vlocN = momyN/rholocalN
         wlocN = momzN/rholocalN

      else

         ulocN = ux
         vlocN = uy
         wlocN = 0.0

      endif


!=======================================================================
! Construct ghost populations
!=======================================================================

      do l = 1,nl


!#######################################################################
! LOWER FACE: j=0
!#######################################################################

         if (j0_is_phys) then

            fint0 = f(l,i,1,k)


!-----------------------------------------------------------------------
! Default: zero-gradient
!-----------------------------------------------------------------------

            f(l,i,0,k) = fint0


!-----------------------------------------------------------------------
! Only cy>0 populations can stream from j=0 into the domain.
!-----------------------------------------------------------------------

            if (cys(l) > 0) then


!=======================================================================
! LOWER INFLOW
!=======================================================================

               if (win0 > 0.0) then

                  m = bounce(l)

                  dfeq = fequil_difference( &
                             l,rholoc0,ux,uy,0.0, &
                             inv1cs2,inv6cs6,ibgk)

                  fin0 = f(m,i,1,k) + dfeq


                  f(l,i,0,k) = &
                       win0*fin0 + &
                       (1.0-win0)*fint0

               endif


!=======================================================================
! LOWER OUTFLOW
!
! NEE:
!
!   fout =
!
!       feq(rho0,u_local)
!     + fint
!     - feq(rho_local,u_local)
!=======================================================================

               if (wout0 > 0.0) then

                  feq_int = fequil_value( &
                               l,rholocal0, &
                               uloc0,vloc0,wloc0, &
                               inv1cs2,inv2cs4,inv6cs6,ibgk)


                  feq_out = fequil_value( &
                               l,rho0, &
                               uloc0,vloc0,wloc0, &
                               inv1cs2,inv2cs4,inv6cs6,ibgk)


                  fout0 = feq_out + &
                          (fint0-feq_int)


                  f(l,i,0,k) = &
                       wout0*fout0 + &
                       (1.0-wout0)*fint0

               endif

            endif

         endif


!#######################################################################
! UPPER FACE: j=ny+1
!#######################################################################

         if (jN_is_phys) then

            fintN = f(l,i,ny,k)


!-----------------------------------------------------------------------
! Default: zero-gradient
!-----------------------------------------------------------------------

            f(l,i,ny+1,k) = fintN


!-----------------------------------------------------------------------
! Only cy<0 populations can stream from j=ny+1 into the domain.
!-----------------------------------------------------------------------

            if (cys(l) < 0) then


!=======================================================================
! UPPER INFLOW
!=======================================================================

               if (winN > 0.0) then

                  m = bounce(l)

                  dfeq = fequil_difference( &
                             l,rholocN,ux,uy,0.0, &
                             inv1cs2,inv6cs6,ibgk)

                  finN = f(m,i,ny,k) + dfeq


                  f(l,i,ny+1,k) = &
                       winN*finN + &
                       (1.0-winN)*fintN

               endif


!=======================================================================
! UPPER OUTFLOW
!=======================================================================

               if (woutN > 0.0) then

                  feq_int = fequil_value( &
                               l,rholocalN, &
                               ulocN,vlocN,wlocN, &
                               inv1cs2,inv2cs4,inv6cs6,ibgk)


                  feq_out = fequil_value( &
                               l,rho0, &
                               ulocN,vlocN,wlocN, &
                               inv1cs2,inv2cs4,inv6cs6,ibgk)


                  foutN = feq_out + &
                          (fintN-feq_int)


                  f(l,i,ny+1,k) = &
                       woutN*foutN + &
                       (1.0-woutN)*fintN

               endif

            endif

         endif

      enddo


#ifndef _CUDA

   enddo
   enddo

!$OMP END PARALLEL DO

#endif


end subroutine boundary_inflow_j_kernel

end module m_boundary_inflow_j_kernel





!WWmodule m_boundary_inflow_j_kernel
!WWcontains
!WW
!WW#ifdef _CUDA
!WW   attributes(global) &
!WW#endif
!WWsubroutine boundary_inflow_j_kernel(f,uvel,udir,rho0,rho_relax, &
!WW                                    inv1cs2,inv2cs4,inv6cs6,ibgk,&
!WW                                    j0_is_phys,jN_is_phys)
!WW
!WW!-----------------------------------------------------------------------
!WW! General inflow / zero-gradient boundary condition in j direction.
!WW!
!WW! For uy > 0:
!WW!
!WW!     j = 0       : prescribed inflow
!WW!     j = ny+1    : zero-gradient outflow
!WW!
!WW! For uy < 0:
!WW!
!WW!     j = 0       : zero-gradient outflow
!WW!     j = ny+1    : prescribed inflow
!WW!
!WW! Around uy = 0 the inflow reconstruction is smoothly switched off.
!WW!
!WW! Only populations actually ENTERING through a face are reconstructed:
!WW!
!WW!     j=0    : cy > 0
!WW!     j=ny+1 : cy < 0
!WW!
!WW! All remaining ghost populations are zero-gradient copies of the
!WW! adjacent interior population.
!WW!
!WW! Incoming populations use
!WW!
!WW!     f_l = f_bounce(l) + feq_l - feq_bounce(l)
!WW!
!WW! with the same equilibrium order as postcoll_kernel.
!WW!
!WW! j0_is_phys / jN_is_phys prevent writes to internal MPI interfaces.
!WW!-----------------------------------------------------------------------
!WW
!WW#ifdef _CUDA
!WW   use cudafor
!WW#endif
!WW
!WW   use mod_dimensions
!WW   use mod_D3Q27setup, only : nl,cys,bounce
!WW   use m_fequil_difference, only : fequil_difference
!WW
!WW   implicit none
!WW
!WW   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
!WW   real, intent(in)    :: uvel(nz)
!WW
!WW   real, value    :: udir
!WW   real, value    :: rho0
!WW   real, value    :: rho_relax
!WW   real, value    :: inv1cs2
!WW   real, value    :: inv6cs6
!WW
!WW   integer, value :: ibgk
!WW
!WW   logical, value :: j0_is_phys
!WW   logical, value :: jN_is_phys
!WW
!WW   integer :: i,k,l,m
!WW
!WW   real, parameter :: pi = acos(-1.0)
!WW   real, parameter :: blend_width = 0.05
!WW
!WW   real :: uxdir,uydir
!WW   real :: uu
!WW   real :: ux,uy
!WW
!WW   real :: rholocal0,rholocalN
!WW   real :: rholoc0,rholocN
!WW
!WW   real :: s
!WW   real :: win0,winN
!WW
!WW   real :: dfeq
!WW   real :: fbc0,fbcN
!WW   real :: fint0,fintN
!WW
!WW
!WW!=======================================================================
!WW! CUDA / CPU indexing
!WW!=======================================================================
!WW
!WW#ifdef _CUDA
!WW
!WW   i = threadIdx%y + (blockIdx%y-1)*blockDim%y
!WW   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
!WW
!WW   if (i < 1 .or. i > nx) return
!WW   if (k < 1 .or. k > nz) return
!WW
!WW#else
!WW
!WW!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                         &
!WW!$OMP& PRIVATE(i,k,l,m,uxdir,uydir,uu,ux,uy,                       &
!WW!$OMP&         rholocal0,rholocalN,rholoc0,rholocN,                &
!WW!$OMP&         s,win0,winN,dfeq,fbc0,fbcN,fint0,fintN)             &
!WW!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,                          &
!WW!$OMP&        inv1cs2,inv6cs6,ibgk,j0_is_phys,jN_is_phys)
!WW
!WW   do k = 1,nz
!WW   do i = 1,nx
!WW
!WW#endif
!WW
!WW
!WW!=======================================================================
!WW! Prescribed velocity
!WW!=======================================================================
!WW
!WW      uxdir = cos(udir*pi/180.0)
!WW      uydir = sin(udir*pi/180.0)
!WW
!WW      uu = uvel(k)
!WW
!WW      ux = uu*uxdir
!WW      uy = uu*uydir
!WW
!WW
!WW!=======================================================================
!WW! Smooth inflow weights
!WW!=======================================================================
!WW
!WW      s = min(1.0,max(0.0,uydir/blend_width))
!WW      win0 = s*s*(3.0-2.0*s)
!WW
!WW      s = min(1.0,max(0.0,-uydir/blend_width))
!WW      winN = s*s*(3.0-2.0*s)
!WW
!WW
!WW!=======================================================================
!WW! Density adjacent to lower j face
!WW!=======================================================================
!WW
!WW      rholocal0 = 0.0
!WW
!WW      do l = 1,nl
!WW         rholocal0 = rholocal0 + f(l,i,1,k)
!WW      enddo
!WW
!WW      rholoc0 = rho_relax*rholocal0 + &
!WW                (1.0-rho_relax)*rho0
!WW
!WW
!WW!=======================================================================
!WW! Density adjacent to upper j face
!WW!=======================================================================
!WW
!WW      rholocalN = 0.0
!WW
!WW      do l = 1,nl
!WW         rholocalN = rholocalN + f(l,i,ny,k)
!WW      enddo
!WW
!WW      rholocN = rho_relax*rholocalN + &
!WW                (1.0-rho_relax)*rho0
!WW
!WW
!WW!=======================================================================
!WW! Construct ghost populations
!WW!=======================================================================
!WW
!WW      do l = 1,nl
!WW
!WW
!WW!-----------------------------------------------------------------------
!WW! LOWER FACE, j=0
!WW!-----------------------------------------------------------------------
!WW
!WW         if (j0_is_phys) then
!WW
!WW            fint0 = f(l,i,1,k)
!WW            fbc0  = fint0
!WW
!WW
!WW!-----------------------------------------------------------------------
!WW! Only cy>0 populations enter through j=0.
!WW!-----------------------------------------------------------------------
!WW
!WW            if (cys(l) > 0 .and. win0 > 0.0) then
!WW
!WW               m = bounce(l)
!WW
!WW               dfeq = fequil_difference( &
!WW                          l,rholoc0,ux,uy,0.0, &
!WW                          inv1cs2,inv6cs6,ibgk)
!WW
!WW               fbc0 = f(m,i,1,k) + dfeq
!WW
!WW            endif
!WW
!WW
!WW            f(l,i,0,k) = &
!WW                 win0*fbc0 + &
!WW                 (1.0-win0)*fint0
!WW
!WW         endif
!WW
!WW
!WW!-----------------------------------------------------------------------
!WW! UPPER FACE, j=ny+1
!WW!-----------------------------------------------------------------------
!WW
!WW         if (jN_is_phys) then
!WW
!WW            fintN = f(l,i,ny,k)
!WW            fbcN  = fintN
!WW
!WW
!WW!-----------------------------------------------------------------------
!WW! Only cy<0 populations enter through j=ny+1.
!WW!-----------------------------------------------------------------------
!WW
!WW            if (cys(l) < 0 .and. winN > 0.0) then
!WW
!WW               m = bounce(l)
!WW
!WW               dfeq = fequil_difference( &
!WW                          l,rholocN,ux,uy,0.0, &
!WW                          inv1cs2,inv6cs6,ibgk)
!WW
!WW               fbcN = f(m,i,ny,k) + dfeq
!WW
!WW            endif
!WW
!WW
!WW            f(l,i,ny+1,k) = &
!WW                 winN*fbcN + &
!WW                 (1.0-winN)*fintN
!WW
!WW         endif
!WW
!WW      enddo
!WW
!WW
!WW#ifndef _CUDA
!WW
!WW   enddo
!WW   enddo
!WW
!WW!$OMP END PARALLEL DO
!WW
!WW#endif
!WW
!WW
!WWend subroutine boundary_inflow_j_kernel
!WW
!WWend module m_boundary_inflow_j_kernel



!!!module m_boundary_inflow_j_kernel
!!!contains
!!!
!!!#ifdef _CUDA
!!!   attributes(global) &
!!!#endif
!!!subroutine boundary_inflow_j_kernel(f,uvel,udir,rho0,rho_relax,j0_is_phys,jN_is_phys)
!!!
!!!!-----------------------------------------------------------------------
!!!! General inflow/outflow boundary condition in the j-direction.
!!!!
!!!! The boundary condition is valid for arbitrary horizontal wind
!!!! direction and provides a smooth transition through uy = 0
!!!! (udir = 0, 180, 360 degrees).
!!!!
!!!! For uy > 0:
!!!!
!!!!     j = 0       : Kruger inflow
!!!!     j = ny+1    : zero-gradient outflow
!!!!
!!!! For uy < 0:
!!!!
!!!!     j = 0       : zero-gradient outflow
!!!!     j = ny+1    : Kruger inflow
!!!!
!!!! Around uy = 0 the Kruger inflow condition is smoothly switched off.
!!!! At exactly uy = 0 both boundaries become zero-normal-gradient:
!!!!
!!!!     f(:,i,0)    = f(:,i,1)
!!!!     f(:,i,ny+1) = f(:,i,ny)
!!!!
!!!! This avoids the previous history-dependent convective outflow
!!!! condition. In particular, no old ghost-cell value is retained when
!!!! a boundary changes from inflow to outflow or vice versa.
!!!!
!!!! j0_is_phys and jN_is_phys determine whether the corresponding
!!!! boundary belongs to the physical domain. This is useful for MPI
!!!! decomposition in the j-direction: internal MPI interfaces are left
!!!! untouched by this routine.
!!!!
!!!! blend_width is expressed in terms of the direction cosine uydir.
!!!! For blend_width = 0.05 the transition is approximately +/- 2.87 deg
!!!! around a direction for which uy = 0.
!!!!
!!!!-----------------------------------------------------------------------
!!!
!!!#ifdef _CUDA
!!!   use cudafor
!!!#endif
!!!
!!!   use mod_dimensions
!!!   use mod_D3Q27setup
!!!
!!!   implicit none
!!!
!!!   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
!!!   real, intent(in)    :: uvel(nz)
!!!
!!!   real, value         :: udir
!!!   real, value         :: rho0
!!!   real, value         :: rho_relax
!!!
!!!   logical, value      :: j0_is_phys
!!!   logical, value      :: jN_is_phys
!!!
!!!   integer :: i,k,l
!!!
!!!   real, parameter :: pi = acos(-1.0)
!!!
!!!   ! Width of transition around uy = 0.
!!!   !
!!!   ! 0.05 corresponds approximately to +/- 2.87 degrees.
!!!   real, parameter :: blend_width = 0.05
!!!
!!!   real :: wl
!!!   real :: cxl,cyl
!!!   real :: uu
!!!
!!!   real :: uxdir,uydir
!!!
!!!   real :: rholocal0,rholocalN
!!!   real :: rholoc0,rholocN
!!!
!!!   real :: fraw0(nl),frawN(nl)
!!!   real :: fkruger0(nl),fkrugerN(nl)
!!!
!!!   real :: s
!!!   real :: win0,winN
!!!
!!!
!!!!=======================================================================
!!!! Thread/grid indexing
!!!!=======================================================================
!!!
!!!#ifdef _CUDA
!!!
!!!   i = threadIdx%y + (blockIdx%y-1)*blockDim%y
!!!   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
!!!
!!!   if (i < 1 .or. i > nx) return
!!!   if (k < 1 .or. k > nz) return
!!!
!!!#else
!!!
!!!!$OMP PARALLEL DO COLLAPSE(2)                                      &
!!!!$OMP& PRIVATE(i,k,l,wl,cxl,cyl,uu,                               &
!!!!$OMP&         uxdir,uydir,                                       &
!!!!$OMP&         rholocal0,rholocalN,rholoc0,rholocN,               &
!!!!$OMP&         fraw0,frawN,fkruger0,fkrugerN,                     &
!!!!$OMP&         s,win0,winN)                                       &
!!!!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,j0_is_phys,jN_is_phys)
!!!
!!!   do k = 1,nz
!!!   do i = 1,nx
!!!
!!!#endif
!!!
!!!
!!!!=======================================================================
!!!! Prescribed horizontal velocity direction
!!!!=======================================================================
!!!
!!!      uxdir = cos(udir*pi/180.0)
!!!      uydir = sin(udir*pi/180.0)
!!!
!!!      uu = uvel(k)
!!!
!!!
!!!!=======================================================================
!!!! Smooth inflow weights
!!!!
!!!! Lower j boundary, j = 0:
!!!!
!!!!        uy <= 0              win0 = 0
!!!!        0 < uy < width       smooth transition
!!!!        uy >= width          win0 = 1
!!!!
!!!! Upper j boundary, j = ny+1:
!!!!
!!!!        uy >= 0              winN = 0
!!!!       -width < uy < 0       smooth transition
!!!!        uy <= -width         winN = 1
!!!!
!!!! Thus at uy = 0:
!!!!
!!!!        win0 = 0
!!!!        winN = 0
!!!!
!!!! and both physical j-boundaries are zero-gradient.
!!!!=======================================================================
!!!
!!!      s = min(1.0,max(0.0,uydir/blend_width))
!!!
!!!      win0 = s*s*(3.0-2.0*s)
!!!
!!!
!!!      s = min(1.0,max(0.0,-uydir/blend_width))
!!!
!!!      winN = s*s*(3.0-2.0*s)
!!!
!!!
!!!!=======================================================================
!!!! Density used in Kruger reconstruction.
!!!!
!!!! The density is evaluated at the first interior plane adjacent to
!!!! each physical boundary and relaxed toward rho0.
!!!!=======================================================================
!!!
!!!      rholocal0 = 0.0
!!!
!!!      do l = 1,nl
!!!         rholocal0 = rholocal0 + f(l,i,1,k)
!!!      enddo
!!!
!!!      rholoc0 = rho_relax*rholocal0 + &
!!!                (1.0-rho_relax)*rho0
!!!
!!!
!!!      rholocalN = 0.0
!!!
!!!      do l = 1,nl
!!!         rholocalN = rholocalN + f(l,i,ny,k)
!!!      enddo
!!!
!!!      rholocN = rho_relax*rholocalN + &
!!!                (1.0-rho_relax)*rho0
!!!
!!!
!!!!=======================================================================
!!!! Raw Kruger reconstruction.
!!!!
!!!! The full prescribed horizontal velocity vector is used:
!!!!
!!!!        u = uu * uxdir
!!!!        v = uu * uydir
!!!!
!!!! This remains important even when the normal velocity component
!!!! through the j-boundary is relatively small.
!!!!=======================================================================
!!!
!!!      do l = 1,nl
!!!
!!!         wl  = weights(l)
!!!
!!!         cxl = real(cxs(l))
!!!         cyl = real(cys(l))
!!!
!!!         fraw0(l) = f(l,i,1,k)                         &
!!!                    - 2.0*wl*rholoc0                  &
!!!                    *(cxl*uu*uxdir + cyl*uu*uydir)    &
!!!                    /cs2
!!!
!!!         frawN(l) = f(l,i,ny,k)                        &
!!!                    - 2.0*wl*rholocN                  &
!!!                    *(cxl*uu*uxdir + cyl*uu*uydir)    &
!!!                    /cs2
!!!
!!!      enddo
!!!
!!!
!!!!=======================================================================
!!!! Construct Kruger ghost populations.
!!!!
!!!! j = 0:
!!!! populations entering the domain have cy > 0.
!!!!
!!!! j = ny+1:
!!!! populations entering the domain have cy < 0.
!!!!
!!!! The bounce mapping reconstructs the incoming populations from their
!!!! opposite velocity directions.
!!!!=======================================================================
!!!
!!!      do l = 1,nl
!!!
!!!         !------------------------------------------------------------
!!!         ! Lower boundary, j = 0
!!!         !------------------------------------------------------------
!!!
!!!         if (cys(l) <= 0) then
!!!            fkruger0(l) = fraw0(l)
!!!         else
!!!            fkruger0(l) = fraw0(bounce(l))
!!!         endif
!!!
!!!
!!!         !------------------------------------------------------------
!!!         ! Upper boundary, j = ny+1
!!!         !------------------------------------------------------------
!!!
!!!         if (cys(l) >= 0) then
!!!            fkrugerN(l) = frawN(l)
!!!         else
!!!            fkrugerN(l) = frawN(bounce(l))
!!!         endif
!!!
!!!      enddo
!!!
!!!
!!!!=======================================================================
!!!! Final ghost-cell boundary condition.
!!!!
!!!! Lower boundary:
!!!!
!!!!   win0 = 1 : Kruger inflow
!!!!   win0 = 0 : zero-normal-gradient
!!!!
!!!! Upper boundary:
!!!!
!!!!   winN = 1 : Kruger inflow
!!!!   winN = 0 : zero-normal-gradient
!!!!
!!!! At uy = 0:
!!!!
!!!!   f(:,i,0)    = f(:,i,1)
!!!!   f(:,i,ny+1) = f(:,i,ny)
!!!!
!!!! Only physical boundaries are modified. Internal MPI boundaries
!!!! remain untouched.
!!!!=======================================================================
!!!
!!!      if (j0_is_phys) then
!!!
!!!         do l = 1,nl
!!!
!!!            f(l,i,0,k) = &
!!!                 win0*fkruger0(l) + &
!!!                 (1.0-win0)*f(l,i,1,k)
!!!
!!!         enddo
!!!
!!!      endif
!!!
!!!
!!!      if (jN_is_phys) then
!!!
!!!         do l = 1,nl
!!!
!!!            f(l,i,ny+1,k) = &
!!!                 winN*fkrugerN(l) + &
!!!                 (1.0-winN)*f(l,i,ny,k)
!!!
!!!         enddo
!!!
!!!      endif
!!!
!!!
!!!#ifndef _CUDA
!!!
!!!   enddo
!!!   enddo
!!!
!!!!$OMP END PARALLEL DO
!!!
!!!#endif
!!!
!!!
!!!end subroutine boundary_inflow_j_kernel
!!!
!!!end module m_boundary_inflow_j_kernel
!!!
!!!


!! module m_boundary_inflow_j_kernel
!! contains
!! #ifdef _CUDA
!!    attributes(global) &
!! #endif
!! subroutine boundary_inflow_j_kernel(f,uvel,udir,rho0,rho_relax,uvel_ref,j0_is_phys,jN_is_phys)
!! ! Inflow/outflow boundary conditions in the j-direction.
!! ! See boundary_inflow_i_kernel for the smoothstep-blend rationale.
!! 
!! #ifdef _CUDA
!!    use cudafor
!! #endif
!!    use mod_dimensions
!!    use mod_D3Q27setup
!! 
!!    implicit none
!!    real, intent(inout) :: f     (nl,0:nx+1,0:ny+1,0:nz+1)
!!    real, intent(in)    :: uvel  (nz)
!!    real, value         :: udir
!!    real, value         :: rho0
!!    real, value         :: rho_relax
!!    real, value         :: uvel_ref
!!    logical, value       :: j0_is_phys
!!    logical, value       :: jN_is_phys
!! 
!!    integer :: i,k,l
!!    real, parameter :: pi = acos(-1.0)
!!    real, parameter :: blend_width    = 0.05 !0.05
!!    real, parameter :: uconv_min_frac = 0.001
!!    real :: wl, cxl, cyl
!!    real :: uu
!!    real :: fraw0(nl), frawN(nl)
!!    real :: fkruger0(nl), fkrugerN(nl)
!!    real :: fconv0(nl), fconvN(nl)
!!    real :: uxdir, uydir
!!    real :: rholoc0, rholocN, rholocal0, rholocalN
!!    real :: uconv, cconv, invden
!!    real :: xi, t, w
!! 
!! !------------------ Indexing (CUDA vs CPU) -------------------------
!! #ifdef _CUDA
!!    i = threadIdx%y + (blockIdx%y-1)*blockDim%y
!!    k = threadIdx%z + (blockIdx%z-1)*blockDim%z
!!    if (i < 1 .or. i > nx) return
!!    if (k < 1 .or. k > nz) return
!! #else
!! !$OMP PARALLEL DO COLLAPSE(2) &
!! !$OMP& PRIVATE(i,k,l,wl,cxl,cyl,uu,fraw0,frawN,fkruger0,fkrugerN,fconv0,fconvN, &
!! !$OMP&         uxdir,uydir,rholoc0,rholocN,rholocal0,rholocalN,uconv,cconv,invden,xi,t,w) &
!! !$OMP& SHARED(f,uvel,udir,rho0,rho_relax,uvel_ref)
!!    do k = 1, nz
!!    do i = 1, nx
!! #endif
!!       uxdir = cos(udir*pi/180.0)
!!       uydir = sin(udir*pi/180.0)
!!       uu    = uvel(k)
!! 
!!       uconv  = max(uu*abs(uydir), uconv_min_frac*uvel_ref)
!!       cconv  = min(1.0, max(0.0, uconv))
!!       invden = 1.0/(1.0+cconv)
!! 
!!       xi = min(1.0, max(-1.0, uydir/blend_width))
!!       t  = 0.5*(xi+1.0)
!!       w  = t*t*(3.0-2.0*t)
!! 
!!       rholocal0 = 0.0
!!       do l = 1, nl
!!          rholocal0 = rholocal0 + f(l,i,1,k)
!!       enddo
!!       rholoc0 = rho_relax*rholocal0 + (1.0-rho_relax)*rho0
!! 
!!       rholocalN = 0.0
!!       do l = 1, nl
!!          rholocalN = rholocalN + f(l,i,ny,k)
!!       enddo
!!       rholocN = rho_relax*rholocalN + (1.0-rho_relax)*rho0
!! 
!!       do l = 1, nl
!!          wl  = weights(l)
!!          cxl = real(cxs(l))
!!          cyl = real(cys(l))
!!          fraw0(l) = f(l,i,1,k)  - 2.0*wl*rholoc0*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
!!          frawN(l) = f(l,i,ny,k) - 2.0*wl*rholocN*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
!!       enddo
!! 
!!       do l = 1, nl
!!          if (cys(l) <= 0) then
!!             fkruger0(l) = fraw0(l)
!!          else
!!             fkruger0(l) = fraw0(bounce(l))
!!          endif
!! 
!!          if (cys(l) >= 0) then
!!             fkrugerN(l) = frawN(l)
!!          else
!!             fkrugerN(l) = frawN(bounce(l))
!!          endif
!!       enddo
!! 
!!       do l = 1, nl
!!          fconv0(l) = (f(l,i,0,k)    + cconv*f(l,i,1,k))  * invden
!!          fconvN(l) = (f(l,i,ny+1,k) + cconv*f(l,i,ny,k)) * invden
!!       enddo
!! 
!!       if (j0_is_phys) then
!!          do l = 1, nl
!!             f(l,i,0,k) = w*fkruger0(l) + (1.0-w)*fconv0(l)
!!          enddo
!!       endif
!!       if (jN_is_phys) then
!!          do l = 1, nl
!!             f(l,i,ny+1,k) = (1.0-w)*fkrugerN(l) + w*fconvN(l)
!!          enddo
!!       endif
!! 
!! #ifndef _CUDA
!!    enddo
!!    enddo
!! !$OMP END PARALLEL DO
!! #endif
!! 
!! end subroutine
!! end module
