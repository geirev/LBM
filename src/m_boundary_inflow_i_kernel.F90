module m_boundary_inflow_i_kernel
contains




!=======================================================================
! i-direction open boundary kernel
!=======================================================================

#ifdef _CUDA
attributes(global) &
#endif
subroutine boundary_inflow_i_kernel(f,uvel,udir,rho0,rho_relax, &
                                    inv1cs2,inv2cs4,inv6cs6,ibgk)

!-----------------------------------------------------------------------
! General inflow / outflow boundary condition in the i direction.
!
!
! INFLOW
! ======
!
! Incoming populations are reconstructed using
!
!       f_l = f_bounce(l) + feq_l - feq_bounce(l)
!
! using the same equilibrium order as postcoll_kernel.
!
!
! OUTFLOW
! =======
!
! Populations entering the computational domain from an OUTLET are
! reconstructed using non-equilibrium extrapolation:
!
!     f_outlet =
!
!         feq(rho0,u_int)
!
!       + [ f_int - feq(rho_int,u_int) ]
!
!
! Thus:
!
!   - rho0 anchors the equilibrium/pressure part at the outlet;
!
!   - the local interior velocity is retained, so a turbine wake can
!     leave the domain without being replaced by the prescribed
!     undisturbed inflow velocity;
!
!   - the local non-equilibrium part is extrapolated through the outlet;
!
!   - only populations streaming FROM the ghost plane INTO the domain
!     require reconstruction.
!
!
! LEFT i FACE:
!
!   cx > 0 populations enter the domain from i=0.
!
!   ux > 0 : inflow reconstruction
!   ux < 0 : NEE outflow reconstruction
!
!
! RIGHT i FACE:
!
!   cx < 0 populations enter the domain from i=nx+1.
!
!   ux < 0 : inflow reconstruction
!   ux > 0 : NEE outflow reconstruction
!
!
! Around ux=0 both inflow and outflow reconstruction are smoothly
! switched off and the boundary approaches zero-gradient.
!
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


!=======================================================================
! Local variables
!=======================================================================

   integer :: j,k,l,m

   real, parameter :: pi = acos(-1.0)

   real, parameter :: blend_width = 0.05


!-----------------------------------------------------------------------
! Prescribed inflow velocity
!-----------------------------------------------------------------------

   real :: uxdir,uydir
   real :: uu
   real :: ux,uy


!-----------------------------------------------------------------------
! Local density at adjacent interior nodes
!-----------------------------------------------------------------------

   real :: rholocal0
   real :: rholocalN

   real :: rholoc0
   real :: rholocN


!-----------------------------------------------------------------------
! Local momentum / velocity at adjacent interior nodes.
!
! Used by the outflow NEE condition.
!-----------------------------------------------------------------------

   real :: momx0,momy0,momz0
   real :: momxN,momyN,momzN

   real :: uloc0,vloc0,wloc0
   real :: ulocN,vlocN,wlocN


!-----------------------------------------------------------------------
! Smooth boundary weights
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

   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z

   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return

#else

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                         &
!$OMP& PRIVATE(j,k,l,m,uxdir,uydir,uu,ux,uy,                       &
!$OMP&         rholocal0,rholocalN,rholoc0,rholocN,                &
!$OMP&         momx0,momy0,momz0,momxN,momyN,momzN,                &
!$OMP&         uloc0,vloc0,wloc0,ulocN,vlocN,wlocN,                &
!$OMP&         s,win0,winN,wout0,woutN,                            &
!$OMP&         fint0,fintN,fin0,finN,fout0,foutN,                  &
!$OMP&         feq_int,feq_out,dfeq)                               &
!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,                          &
!$OMP&        inv1cs2,inv2cs4,inv6cs6,ibgk)

   do k = 1,nz
   do j = 1,ny

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
! Smooth inflow/outflow activation weights.
!
! LEFT:
!
!     win0  -> 1 for uxdir > +blend_width
!     wout0 -> 1 for uxdir < -blend_width
!
! RIGHT:
!
!     winN  -> 1 for uxdir < -blend_width
!     woutN -> 1 for uxdir > +blend_width
!
! At uxdir=0 all four weights are zero.
!=======================================================================

      s = min(1.0,max(0.0, uxdir/blend_width))
      win0  = s*s*(3.0-2.0*s)
      woutN = win0


      s = min(1.0,max(0.0,-uxdir/blend_width))
      winN  = s*s*(3.0-2.0*s)
      wout0 = winN


!=======================================================================
! Density and momentum at adjacent interior nodes.
!
! These moments are evaluated directly from f because the boundary is
! applied to the current distribution state.
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

         rholocal0 = rholocal0 + f(l,1,j,k)

         momx0 = momx0 + real(cxs(l))*f(l,1,j,k)
         momy0 = momy0 + real(cys(l))*f(l,1,j,k)
         momz0 = momz0 + real(czs(l))*f(l,1,j,k)


         rholocalN = rholocalN + f(l,nx,j,k)

         momxN = momxN + real(cxs(l))*f(l,nx,j,k)
         momyN = momyN + real(cys(l))*f(l,nx,j,k)
         momzN = momzN + real(czs(l))*f(l,nx,j,k)

      enddo


!=======================================================================
! Density used by prescribed inflow reconstruction.
!=======================================================================

      rholoc0 = rho_relax*rholocal0 + &
                (1.0-rho_relax)*rho0

      rholocN = rho_relax*rholocalN + &
                (1.0-rho_relax)*rho0


!=======================================================================
! Local velocity used by NEE outflow reconstruction.
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
! LEFT FACE: i=0
!#######################################################################

         fint0 = f(l,1,j,k)


!-----------------------------------------------------------------------
! Default: zero-gradient.
!-----------------------------------------------------------------------

         f(l,0,j,k) = fint0


!-----------------------------------------------------------------------
! Only cx>0 populations can stream from i=0 back into the domain.
!-----------------------------------------------------------------------

         if (cxs(l) > 0) then


!=======================================================================
! LEFT INFLOW
!=======================================================================

            if (win0 > 0.0) then

               m = bounce(l)

               dfeq = fequil_difference( &
                          l,rholoc0,ux,uy,0.0, &
                          inv1cs2,inv6cs6,ibgk)

               fin0 = f(m,1,j,k) + dfeq


               f(l,0,j,k) = &
                    win0*fin0 + &
                    (1.0-win0)*fint0

            endif


!=======================================================================
! LEFT OUTFLOW
!
! Non-equilibrium extrapolation:
!
! fout =
!
!   feq(rho0,u_local)
! + fint
! - feq(rho_local,u_local)
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


               f(l,0,j,k) = &
                    wout0*fout0 + &
                    (1.0-wout0)*fint0

            endif

         endif


!#######################################################################
! RIGHT FACE: i=nx+1
!#######################################################################

         fintN = f(l,nx,j,k)


!-----------------------------------------------------------------------
! Default: zero-gradient.
!-----------------------------------------------------------------------

         f(l,nx+1,j,k) = fintN


!-----------------------------------------------------------------------
! Only cx<0 populations can stream from i=nx+1 back into the domain.
!-----------------------------------------------------------------------

         if (cxs(l) < 0) then


!=======================================================================
! RIGHT INFLOW
!=======================================================================

            if (winN > 0.0) then

               m = bounce(l)

               dfeq = fequil_difference( &
                          l,rholocN,ux,uy,0.0, &
                          inv1cs2,inv6cs6,ibgk)

               finN = f(m,nx,j,k) + dfeq


               f(l,nx+1,j,k) = &
                    winN*finN + &
                    (1.0-winN)*fintN

            endif


!=======================================================================
! RIGHT OUTFLOW
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


               f(l,nx+1,j,k) = &
                    woutN*foutN + &
                    (1.0-woutN)*fintN

            endif

         endif

      enddo


#ifndef _CUDA

   enddo
   enddo

!$OMP END PARALLEL DO

#endif


end subroutine boundary_inflow_i_kernel

end module m_boundary_inflow_i_kernel


!WW! module m_boundary_inflow_i_kernel
!WW! contains
!WW! 
!WW! #ifdef _CUDA
!WW!    attributes(global) &
!WW! #endif
!WW! subroutine boundary_inflow_i_kernel(f,uvel,udir,rho0,rho_relax, &
!WW!                                     inv1cs2,inv6cs6,ibgk)
!WW! 
!WW! !-----------------------------------------------------------------------
!WW! ! General inflow / zero-gradient boundary condition in i direction.
!WW! !
!WW! ! For ux > 0:
!WW! !
!WW! !     i = 0       : prescribed inflow
!WW! !     i = nx+1    : zero-gradient outflow
!WW! !
!WW! ! For ux < 0:
!WW! !
!WW! !     i = 0       : zero-gradient outflow
!WW! !     i = nx+1    : prescribed inflow
!WW! !
!WW! ! Around ux = 0 the inflow reconstruction is smoothly switched off.
!WW! !
!WW! ! IMPORTANT:
!WW! !
!WW! ! Only populations that actually ENTER the computational domain through
!WW! ! a given face are reconstructed.
!WW! !
!WW! !     i=0    : cx > 0
!WW! !     i=nx+1 : cx < 0
!WW! !
!WW! ! All other populations are copied from the adjacent interior node.
!WW! !
!WW! ! The incoming population is reconstructed using
!WW! !
!WW! !     f_l = f_bounce(l) + feq_l - feq_bounce(l)
!WW! !
!WW! ! where feq_l-feq_bounce(l) is evaluated consistently with the
!WW! ! equilibrium used in postcoll_kernel, including the third-order
!WW! ! Hermite contribution when ibgk == 3.
!WW! !
!WW! ! Consequently an exact uniform equilibrium is a fixed point of this
!WW! ! boundary treatment.
!WW! !-----------------------------------------------------------------------
!WW! 
!WW! #ifdef _CUDA
!WW!    use cudafor
!WW! #endif
!WW! 
!WW!    use mod_dimensions
!WW!    use mod_D3Q27setup, only : nl,cxs,bounce
!WW!    use m_fequil_difference, only : fequil_difference
!WW! 
!WW!    implicit none
!WW! 
!WW!    real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
!WW!    real, intent(in)    :: uvel(nz)
!WW! 
!WW!    real, value    :: udir
!WW!    real, value    :: rho0
!WW!    real, value    :: rho_relax
!WW!    real, value    :: inv1cs2
!WW!    real, value    :: inv6cs6
!WW! 
!WW!    integer, value :: ibgk
!WW! 
!WW!    integer :: j,k,l,m
!WW! 
!WW!    real, parameter :: pi = acos(-1.0)
!WW!    real, parameter :: blend_width = 0.05
!WW! 
!WW!    real :: uxdir,uydir
!WW!    real :: uu
!WW!    real :: ux,uy
!WW! 
!WW!    real :: rholocal0,rholocalN
!WW!    real :: rholoc0,rholocN
!WW! 
!WW!    real :: s
!WW!    real :: win0,winN
!WW! 
!WW!    real :: dfeq
!WW!    real :: fbc0,fbcN
!WW!    real :: fint0,fintN
!WW! 
!WW! 
!WW! !=======================================================================
!WW! ! CUDA / CPU indexing
!WW! !=======================================================================
!WW! 
!WW! #ifdef _CUDA
!WW! 
!WW!    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
!WW!    k = threadIdx%z + (blockIdx%z-1)*blockDim%z
!WW! 
!WW!    if (j < 1 .or. j > ny) return
!WW!    if (k < 1 .or. k > nz) return
!WW! 
!WW! #else
!WW! 
!WW! !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                         &
!WW! !$OMP& PRIVATE(j,k,l,m,uxdir,uydir,uu,ux,uy,                       &
!WW! !$OMP&         rholocal0,rholocalN,rholoc0,rholocN,                &
!WW! !$OMP&         s,win0,winN,dfeq,fbc0,fbcN,fint0,fintN)             &
!WW! !$OMP& SHARED(f,uvel,udir,rho0,rho_relax,                          &
!WW! !$OMP&        inv1cs2,inv6cs6,ibgk)
!WW! 
!WW!    do k = 1,nz
!WW!    do j = 1,ny
!WW! 
!WW! #endif
!WW! 
!WW! 
!WW! !=======================================================================
!WW! ! Prescribed velocity
!WW! !=======================================================================
!WW! 
!WW!       uxdir = cos(udir*pi/180.0)
!WW!       uydir = sin(udir*pi/180.0)
!WW! 
!WW!       uu = uvel(k)
!WW! 
!WW!       ux = uu*uxdir
!WW!       uy = uu*uydir
!WW! 
!WW! 
!WW! !=======================================================================
!WW! ! Smooth inflow weights
!WW! !
!WW! ! i=0:
!WW! !
!WW! !    ux <= 0            win0 = 0
!WW! !    ux >= blend_width  win0 = 1
!WW! !
!WW! ! i=nx+1:
!WW! !
!WW! !    ux >= 0            winN = 0
!WW! !    ux <= -blend_width winN = 1
!WW! !=======================================================================
!WW! 
!WW!       s = min(1.0,max(0.0,uxdir/blend_width))
!WW!       win0 = s*s*(3.0-2.0*s)
!WW! 
!WW!       s = min(1.0,max(0.0,-uxdir/blend_width))
!WW!       winN = s*s*(3.0-2.0*s)
!WW! 
!WW! 
!WW! !=======================================================================
!WW! ! Density used for equilibrium reconstruction
!WW! !=======================================================================
!WW! 
!WW!       rholocal0 = 0.0
!WW! 
!WW!       do l = 1,nl
!WW!          rholocal0 = rholocal0 + f(l,1,j,k)
!WW!       enddo
!WW! 
!WW!       rholoc0 = rho_relax*rholocal0 + &
!WW!                 (1.0-rho_relax)*rho0
!WW! 
!WW! 
!WW!       rholocalN = 0.0
!WW! 
!WW!       do l = 1,nl
!WW!          rholocalN = rholocalN + f(l,nx,j,k)
!WW!       enddo
!WW! 
!WW!       rholocN = rho_relax*rholocalN + &
!WW!                 (1.0-rho_relax)*rho0
!WW! 
!WW! 
!WW! !=======================================================================
!WW! ! Construct ghost populations
!WW! !=======================================================================
!WW! 
!WW!       do l = 1,nl
!WW! 
!WW! !-----------------------------------------------------------------------
!WW! ! LEFT FACE, i=0
!WW! !
!WW! ! Start with zero-gradient continuation for ALL populations.
!WW! !-----------------------------------------------------------------------
!WW! 
!WW!          fint0 = f(l,1,j,k)
!WW!          fbc0  = fint0
!WW! 
!WW! 
!WW! !-----------------------------------------------------------------------
!WW! ! Only cx>0 populations enter the domain through i=0.
!WW! !-----------------------------------------------------------------------
!WW! 
!WW!          if (cxs(l) > 0 .and. win0 > 0.0) then
!WW! 
!WW!             m = bounce(l)
!WW! 
!WW!             dfeq = fequil_difference( &
!WW!                        l,rholoc0,ux,uy,0.0, &
!WW!                        inv1cs2,inv6cs6,ibgk)
!WW! 
!WW!             fbc0 = f(m,1,j,k) + dfeq
!WW! 
!WW!          endif
!WW! 
!WW! 
!WW! !-----------------------------------------------------------------------
!WW! ! Smooth transition between inflow reconstruction and zero-gradient.
!WW! !-----------------------------------------------------------------------
!WW! 
!WW!          f(l,0,j,k) = &
!WW!               win0*fbc0 + &
!WW!               (1.0-win0)*fint0
!WW! 
!WW! 
!WW! !-----------------------------------------------------------------------
!WW! ! RIGHT FACE, i=nx+1
!WW! !
!WW! ! Start with zero-gradient continuation for ALL populations.
!WW! !-----------------------------------------------------------------------
!WW! 
!WW!          fintN = f(l,nx,j,k)
!WW!          fbcN  = fintN
!WW! 
!WW! 
!WW! !-----------------------------------------------------------------------
!WW! ! Only cx<0 populations enter the domain through i=nx+1.
!WW! !-----------------------------------------------------------------------
!WW! 
!WW!          if (cxs(l) < 0 .and. winN > 0.0) then
!WW! 
!WW!             m = bounce(l)
!WW! 
!WW!             dfeq = fequil_difference( &
!WW!                        l,rholocN,ux,uy,0.0, &
!WW!                        inv1cs2,inv6cs6,ibgk)
!WW! 
!WW!             fbcN = f(m,nx,j,k) + dfeq
!WW! 
!WW!          endif
!WW! 
!WW! 
!WW! !-----------------------------------------------------------------------
!WW! ! Smooth transition between inflow reconstruction and zero-gradient.
!WW! !-----------------------------------------------------------------------
!WW! 
!WW!          f(l,nx+1,j,k) = &
!WW!               winN*fbcN + &
!WW!               (1.0-winN)*fintN
!WW! 
!WW!       enddo
!WW! 
!WW! 
!WW! #ifndef _CUDA
!WW! 
!WW!    enddo
!WW!    enddo
!WW! 
!WW! !$OMP END PARALLEL DO
!WW! 
!WW! #endif
!WW! 
!WW! 
!WW! end subroutine boundary_inflow_i_kernel
!WW! 
!WW! end module m_boundary_inflow_i_kernel
!WW! 




!!!module m_boundary_inflow_i_kernel
!!!contains
!!!
!!!#ifdef _CUDA
!!!   attributes(global) &
!!!#endif
!!!subroutine boundary_inflow_i_kernel(f,uvel,udir,rho0,rho_relax)
!!!
!!!!-----------------------------------------------------------------------
!!!! General inflow/outflow boundary condition in the i-direction.
!!!!
!!!! The boundary condition is designed for arbitrary horizontal wind
!!!! direction and, in particular, for a smooth change of flow direction
!!!! through ux = 0 (udir = 90 or 270 degrees).
!!!!
!!!! For ux > 0:
!!!!
!!!!     i = 0       : Kruger inflow
!!!!     i = nx+1    : zero-gradient outflow
!!!!
!!!! For ux < 0:
!!!!
!!!!     i = 0       : zero-gradient outflow
!!!!     i = nx+1    : Kruger inflow
!!!!
!!!! Around ux = 0 the Kruger inflow condition is smoothly switched off.
!!!! At exactly ux = 0 both boundaries become zero-normal-gradient:
!!!!
!!!!     f(:,0)    = f(:,1)
!!!!     f(:,nx+1) = f(:,nx)
!!!!
!!!! This avoids blending an inflow condition with a history-dependent
!!!! convective outflow condition when the flow direction changes through
!!!! tangential to the boundary.
!!!!
!!!! The inflow weights are independent:
!!!!
!!!!     win0 = 1       full inflow at i=0
!!!!     win0 = 0       zero-gradient at i=0
!!!!
!!!!     winN = 1       full inflow at i=nx+1
!!!!     winN = 0       zero-gradient at i=nx+1
!!!!
!!!! blend_width is expressed in terms of the direction cosine uxdir.
!!!! For blend_width = 0.05 the transition around 90 degrees is
!!!! approximately +/- 2.87 degrees.
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
!!!   integer :: j,k,l
!!!
!!!   real, parameter :: pi = acos(-1.0)
!!!
!!!   ! Width of transition around ux = 0.
!!!   !
!!!   ! 0.05 corresponds approximately to +/- 2.87 degrees
!!!   ! around 90 and 270 degrees.
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
!!!   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
!!!   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
!!!
!!!   if (j < 1 .or. j > ny) return
!!!   if (k < 1 .or. k > nz) return
!!!
!!!#else
!!!
!!!!$OMP PARALLEL DO COLLAPSE(2)                                      &
!!!!$OMP& PRIVATE(j,k,l,wl,cxl,cyl,uu,                               &
!!!!$OMP&         uxdir,uydir,                                       &
!!!!$OMP&         rholocal0,rholocalN,rholoc0,rholocN,               &
!!!!$OMP&         fraw0,frawN,fkruger0,fkrugerN,                     &
!!!!$OMP&         s,win0,winN)                                       &
!!!!$OMP& SHARED(f,uvel,udir,rho0,rho_relax)
!!!
!!!   do k = 1,nz
!!!   do j = 1,ny
!!!
!!!#endif
!!!
!!!
!!!!=======================================================================
!!!! Prescribed velocity
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
!!!! Left boundary:
!!!!
!!!!        ux <= 0              win0 = 0
!!!!        0 < ux < width       smooth transition
!!!!        ux >= width          win0 = 1
!!!!
!!!! Right boundary:
!!!!
!!!!        ux >= 0              winN = 0
!!!!       -width < ux < 0       smooth transition
!!!!        ux <= -width         winN = 1
!!!!
!!!! Thus at ux = 0:
!!!!
!!!!        win0 = 0
!!!!        winN = 0
!!!!
!!!! and both i-boundaries are zero-gradient.
!!!!=======================================================================
!!!
!!!      s = min(1.0,max(0.0,uxdir/blend_width))
!!!
!!!      win0 = s*s*(3.0-2.0*s)
!!!
!!!
!!!      s = min(1.0,max(0.0,-uxdir/blend_width))
!!!
!!!      winN = s*s*(3.0-2.0*s)
!!!
!!!
!!!!=======================================================================
!!!! Density used in Kruger reconstruction
!!!!
!!!! The density is obtained from the first interior plane and optionally
!!!! relaxed toward rho0.
!!!!=======================================================================
!!!
!!!      rholocal0 = 0.0
!!!
!!!      do l = 1,nl
!!!         rholocal0 = rholocal0 + f(l,1,j,k)
!!!      enddo
!!!
!!!      rholoc0 = rho_relax*rholocal0 + &
!!!                (1.0-rho_relax)*rho0
!!!
!!!
!!!      rholocalN = 0.0
!!!
!!!      do l = 1,nl
!!!         rholocalN = rholocalN + f(l,nx,j,k)
!!!      enddo
!!!
!!!      rholocN = rho_relax*rholocalN + &
!!!                (1.0-rho_relax)*rho0
!!!
!!!
!!!!=======================================================================
!!!! Raw Kruger reconstruction
!!!!
!!!! Use the complete prescribed horizontal velocity vector:
!!!!
!!!!        u = uu * uxdir
!!!!        v = uu * uydir
!!!!
!!!! rather than only its normal component.
!!!!=======================================================================
!!!
!!!      do l = 1,nl
!!!
!!!         wl  = weights(l)
!!!
!!!         cxl = real(cxs(l))
!!!         cyl = real(cys(l))
!!!
!!!         fraw0(l) = f(l,1,j,k)                         &
!!!                    - 2.0*wl*rholoc0                  &
!!!                    *(cxl*uu*uxdir + cyl*uu*uydir)    &
!!!                    /cs2
!!!
!!!         frawN(l) = f(l,nx,j,k)                        &
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
!!!! i = 0:
!!!! populations entering the domain have cx > 0.
!!!!
!!!! i = nx+1:
!!!! populations entering the domain have cx < 0.
!!!!
!!!! The bounce mapping reconstructs these populations from their
!!!! opposite directions.
!!!!=======================================================================
!!!
!!!      do l = 1,nl
!!!
!!!         !------------------------------------------------------------
!!!         ! Left boundary, i = 0
!!!         !------------------------------------------------------------
!!!
!!!
!!!         if (cxs(l) <= 0) then
!!!            fkruger0(l) = fraw0(l)
!!!         else
!!!            fkruger0(l) = fraw0(bounce(l))
!!!         endif
!!!
!!!
!!!         !------------------------------------------------------------
!!!         ! Right boundary, i = nx+1
!!!         !------------------------------------------------------------
!!!
!!!         if (cxs(l) >= 0) then
!!!
!!!            fkrugerN(l) = frawN(l)
!!!
!!!         else
!!!
!!!            fkrugerN(l) = frawN(bounce(l))
!!!
!!!         endif
!!!
!!!      enddo
!!!
!!!
!!!!=======================================================================
!!!! Final ghost-cell boundary condition
!!!!
!!!! Left boundary:
!!!!
!!!!   win0 = 1 : Kruger inflow
!!!!   win0 = 0 : zero-normal-gradient
!!!!
!!!! Right boundary:
!!!!
!!!!   winN = 1 : Kruger inflow
!!!!   winN = 0 : zero-normal-gradient
!!!!
!!!! In particular, at uxdir = 0:
!!!!
!!!!   f(:,0)    = f(:,1)
!!!!   f(:,nx+1) = f(:,nx)
!!!!
!!!! No previous ghost-cell value enters the calculation, so there is no
!!!! memory of the previous inflow/outflow role of the boundary.
!!!!=======================================================================
!!!
!!!      do l = 1,nl
!!!
!!!         f(l,0,j,k) = &
!!!              win0*fkruger0(l) + &
!!!              (1.0-win0)*f(l,1,j,k)
!!!
!!!
!!!         f(l,nx+1,j,k) = &
!!!              winN*fkrugerN(l) + &
!!!              (1.0-winN)*f(l,nx,j,k)
!!!
!!!      enddo
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
!!!end subroutine boundary_inflow_i_kernel
!!!
!!!end module m_boundary_inflow_i_kernel


!! module m_boundary_inflow_i_kernel
!! contains
!! #ifdef _CUDA
!!    attributes(global) &
!! #endif
!! subroutine boundary_inflow_i_kernel(f,uvel,udir,rho0,rho_relax,uvel_ref)
!! ! Inflow/outflow boundary conditions in the i-direction.
!! !
!! !  - Both ghost planes (i=0, i=nx+1) are built every step as a smooth,
!! !    C1-continuous blend of a Kruger-type inflow reconstruction and a
!! !    convective (Orlanski-type) outflow update, weighted by a
!! !    smoothstep function of uxdir. This replaces a hard switch between
!! !    formulas at uxdir=0, which previously injected a discontinuity
!! !    into the ghost layer every time the flow direction rotated
!! !    through tangential to this face - visible downstream as a
!! !    perturbation seeding wake meandering.
!! !  - blend_width sets the angular half-width (in direction-cosine
!! !    units) of the transition zone straddling tangential flow. Wider
!! !    = smoother but blends inflow/outflow physics over a broader
!! !    angular range; narrower = sharper but less smoothing. Tune
!! !    against your rotation rate - it should be wide enough, in
!! !    timesteps, for pressure/momentum to equilibrate as the direction
!! !    sweeps through it.
!! !  - uconv is floored against uvel_ref (a hoisted reference speed,
!! !    e.g. maxval(uvel) computed once outside this kernel) so a small
!! !    or zero local uxdir projection never freezes the outflow ghost.
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
!!    real, value         :: uvel_ref   ! hoisted bulk reference speed, e.g. maxval(uvel)
!! 
!!    integer :: j,k,l
!!    real, parameter :: pi = acos(-1.0)
!!    real, parameter :: blend_width     = 0.05 !0.05   ! tune; direction-cosine half-width
!!    real, parameter :: uconv_min_frac  = 0.001
!!    real :: wl, cxl, cyl
!!    real :: uu
!!    real :: fraw0(nl), frawN(nl)          ! raw Kruger correction at each plane
!!    real :: fkruger0(nl), fkrugerN(nl)    ! bounce-mapped Kruger ghost at each plane
!!    real :: fconv0(nl), fconvN(nl)        ! convective-outflow candidate at each plane
!!    real :: uxdir, uydir
!!    real :: rholoc0, rholocN, rholocal0, rholocalN
!!    real :: uconv, cconv, invden
!!    real :: xi, t, w
!! 
!! 
!! !------------------ Indexing (CUDA vs CPU) -------------------------
!! #ifdef _CUDA
!!    j = threadIdx%y + (blockIdx%y-1)*blockDim%y
!!    k = threadIdx%z + (blockIdx%z-1)*blockDim%z
!!    if (j < 1 .or. j > ny) return
!!    if (k < 1 .or. k > nz) return
!! #else
!! !$OMP PARALLEL DO COLLAPSE(2) &
!! !$OMP& PRIVATE(j,k,l,wl,cxl,cyl,uu,fraw0,frawN,fkruger0,fkrugerN,fconv0,fconvN, &
!! !$OMP&         uxdir,uydir,rholoc0,rholocN,rholocal0,rholocalN,uconv,cconv,invden,xi,t,w) &
!! !$OMP& SHARED(f,uvel,udir,rho0,rho_relax,uvel_ref)
!!    do k = 1, nz
!!    do j = 1, ny
!! #endif
!!       uxdir = cos(udir*pi/180.0)
!!       uydir = sin(udir*pi/180.0)
!!       uu    = uvel(k)
!! 
!!       ! Convective speed: floored so a small/zero i-normal component
!!       ! never freezes whichever plane is currently outflow-weighted.
!!       uconv  = max(uu*abs(uxdir), uconv_min_frac*uvel_ref)
!!       cconv  = min(1.0, max(0.0, uconv))
!!       invden = 1.0/(1.0+cconv)
!! 
!!       ! Smoothstep blend weight: w->1 as uxdir>>0 (i=0 fully inflow,
!!       ! nx+1 fully outflow), w->0 as uxdir<<0 (roles swapped), and a
!!       ! continuous, zero-derivative-at-both-ends transition in between.
!!       xi = min(1.0, max(-1.0, uxdir/blend_width))
!!       t  = 0.5*(xi+1.0)
!!       w  = t*t*(3.0-2.0*t)
!! 
!!       !-----------------------------------------------------------
!!       ! Local density at each plane (relaxed blend, as before).
!!       !-----------------------------------------------------------
!!       rholocal0 = 0.0
!!       do l = 1, nl
!!          rholocal0 = rholocal0 + f(l,1,j,k)
!!       enddo
!!       rholoc0 = rho_relax*rholocal0 + (1.0-rho_relax)*rho0
!! 
!!       rholocalN = 0.0
!!       do l = 1, nl
!!          rholocalN = rholocalN + f(l,nx,j,k)
!!       enddo
!!       rholocN = rho_relax*rholocalN + (1.0-rho_relax)*rho0
!! 
!!       !-----------------------------------------------------------
!!       ! Raw Kruger correction at each plane (direction-general;
!!       ! valid for either sign of uxdir).
!!       !-----------------------------------------------------------
!!       do l = 1, nl
!!          wl  = weights(l)
!!          cxl = real(cxs(l))
!!          cyl = real(cys(l))
!!          fraw0(l) = f(l,1,j,k)  - 2.0*wl*rholoc0*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
!!          frawN(l) = f(l,nx,j,k) - 2.0*wl*rholocN*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
!!       enddo
!! 
!!       !-----------------------------------------------------------
!!       ! Bounce-mapped Kruger ghost set at each plane (same mapping
!!       ! logic as before, now computed unconditionally at both planes).
!!       !-----------------------------------------------------------
!!       do l = 1, nl
!!          if (cxs(l) <= 0) then
!!             fkruger0(l) = fraw0(l)
!!          else
!!             fkruger0(l) = fraw0(bounce(l))
!!          endif
!! 
!!          if (cxs(l) >= 0) then
!!             fkrugerN(l) = frawN(l)
!!          else
!!             fkrugerN(l) = frawN(bounce(l))
!!          endif
!!       enddo
!! 
!!       !-----------------------------------------------------------
!!       ! Convective-outflow candidate at each plane (old ghost value
!!       ! + freshly streamed interior neighbor, as before).
!!       !-----------------------------------------------------------
!!       do l = 1, nl
!!          fconv0(l) = (f(l,0,j,k)    + cconv*f(l,1,j,k))  * invden
!!          fconvN(l) = (f(l,nx+1,j,k) + cconv*f(l,nx,j,k)) * invden
!!       enddo
!! 
!!       !-----------------------------------------------------------
!!       ! Smooth blend at each plane. w=1: i=0 pure inflow, nx+1 pure
!!       ! outflow. w=0: roles swapped. Continuous through w=0.5.
!!       !-----------------------------------------------------------
!!       do l = 1, nl
!!          f(l,0,j,k)    = w*fkruger0(l)          + (1.0-w)*fconv0(l)
!!          f(l,nx+1,j,k) = (1.0-w)*fkrugerN(l)     + w*fconvN(l)
!!       enddo
!! 
!! #ifndef _CUDA
!!    enddo
!!    enddo
!! !$OMP END PARALLEL DO
!! #endif
!! 
!! end subroutine
!! end module
