module m_boundary_inflow_edges_ij_kernel
contains

#ifdef _CUDA
attributes(global) &
#endif
subroutine boundary_inflow_edges_ij_kernel( &
   f,uvel,udir,rho0,rho_relax,inv1cs2,inv6cs6,ibgk, &
   j0_is_phys,jN_is_phys)

!-----------------------------------------------------------------------
! General inflow/zero-gradient treatment of the four true i-j edges
! (intersection of two side boundaries), for k = 1..nz.
!
! This routine must be called AFTER:
!
!    boundary_inflow_i_kernel
!    boundary_inflow_j_kernel
!
! because the baseline edge values use the already completed face
! ghost distributions.
!
!
! BASELINE EDGE VALUE
! -------------------
!
! Each edge is first constructed from the two adjacent face ghost
! distributions:
!
!                  |ux| f_i-face + |uy| f_j-face
!       fbase = ----------------------------------
!                         |ux| + |uy|
!
! Therefore:
!
!    ux -> 0 : edge becomes continuation of the j-face condition
!    uy -> 0 : edge becomes continuation of the i-face condition
!
!
! TRUE INFLOW-INFLOW POPULATIONS
! ------------------------------
!
! Only populations crossing BOTH boundaries are reconstructed:
!
!    (0,0)       : cx > 0, cy > 0
!    (0,ny+1)    : cx > 0, cy < 0
!    (nx+1,0)    : cx < 0, cy > 0
!    (nx+1,ny+1) : cx < 0, cy < 0
!
! These populations use
!
!       f_l = f_bounce(l) + feq(l) - feq(bounce(l))
!
! where feq(l)-feq(bounce(l)) is evaluated using exactly the same
! equilibrium order as postcoll_kernel.
!
! This includes the third-order Hermite contribution when ibgk == 3.
!
!
! SMOOTH ROLE TRANSITION
! ----------------------
!
! The true inflow-inflow edge reconstruction is weighted by the product
! of the corresponding two face inflow weights:
!
!    win00 = win_i0 * win_j0
!    win0N = win_i0 * win_jN
!    winN0 = win_iN * win_j0
!    winNN = win_iN * win_jN
!
! Hence the special edge reconstruction vanishes smoothly if either
! face becomes tangential or outflow.
!
!
! MPI
! ---
!
! j0_is_phys / jN_is_phys prevent modification of edges lying on
! internal MPI j interfaces.
!
!-----------------------------------------------------------------------

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions
   use mod_D3Q27setup, only : nl,cxs,cys,bounce
   use m_fequil_difference, only : fequil_difference

   implicit none

   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: uvel(nz)

   real, value    :: udir
   real, value    :: rho0
   real, value    :: rho_relax
   real, value    :: inv1cs2
   real, value    :: inv6cs6

   integer, value :: ibgk

   logical, value :: j0_is_phys
   logical, value :: jN_is_phys

#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uvel
#endif

   integer :: k,l,m,mm

   real, parameter :: pi          = acos(-1.0)
   real, parameter :: blend_width = 0.05

   real :: uxdir,uydir
   real :: ax,ay,aden

   real :: uu
   real :: ux,uy

   real :: s

   real :: win_i0,win_iN
   real :: win_j0,win_jN

   real :: win00,win0N
   real :: winN0,winNN

   real :: fbase
   real :: fedge

   real :: rholocal
   real :: rhocorner

   real :: dfeq


!=======================================================================
! CUDA / CPU indexing
!=======================================================================

#ifdef _CUDA

   l = threadIdx%x + (blockIdx%x-1)*blockDim%x
   k = threadIdx%y + (blockIdx%y-1)*blockDim%y

   if (l < 1 .or. l > nl) return
   if (k < 1 .or. k > nz) return

#else

!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                          &
!$OMP& PRIVATE(l,k,m,mm,uxdir,uydir,ax,ay,aden,uu,ux,uy,s,          &
!$OMP&         win_i0,win_iN,win_j0,win_jN,                         &
!$OMP&         win00,win0N,winN0,winNN,                             &
!$OMP&         fbase,fedge,rholocal,rhocorner,dfeq)                 &
!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,                           &
!$OMP&        inv1cs2,inv6cs6,ibgk,j0_is_phys,jN_is_phys)

   do k = 1,nz
   do l = 1,nl

#endif


!=======================================================================
! Prescribed velocity
!=======================================================================

      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)

      uu = uvel(k)

      ux = uu*uxdir
      uy = uu*uydir

      ax = abs(uxdir)
      ay = abs(uydir)

      ! uxdir^2 + uydir^2 = 1, so aden is always positive.
      aden = ax + ay


!=======================================================================
! Smooth per-face inflow weights.
!
! Must use exactly the same blend_width and smoothstep as the i and j
! face routines.
!=======================================================================

      s = min(1.0,max(0.0, uxdir/blend_width))
      win_i0 = s*s*(3.0-2.0*s)

      s = min(1.0,max(0.0,-uxdir/blend_width))
      win_iN = s*s*(3.0-2.0*s)

      s = min(1.0,max(0.0, uydir/blend_width))
      win_j0 = s*s*(3.0-2.0*s)

      s = min(1.0,max(0.0,-uydir/blend_width))
      win_jN = s*s*(3.0-2.0*s)


!=======================================================================
! True two-face inflow weights.
!=======================================================================

      win00 = win_i0*win_j0
      win0N = win_i0*win_jN

      winN0 = win_iN*win_j0
      winNN = win_iN*win_jN


!=======================================================================
! Edge (i=0,j=0)
!
! Adjacent completed face values:
!
!    i-face : f(l,0,1,k)
!    j-face : f(l,1,0,k)
!=======================================================================

      if (j0_is_phys) then

         fbase = ( ax*f(l,0,1,k) + &
                   ay*f(l,1,0,k) ) / aden

         f(l,0,0,k) = fbase


!-----------------------------------------------------------------------
! Only populations entering through BOTH i=0 and j=0 are reconstructed:
!
!       cx > 0 and cy > 0
!-----------------------------------------------------------------------

         if (win00 > 0.0 .and. &
             cxs(l) > 0 .and. cys(l) > 0) then

            rholocal = 0.0

            do mm = 1,nl
               rholocal = rholocal + f(mm,1,1,k)
            enddo

            rhocorner = rho_relax*rholocal + &
                        (1.0-rho_relax)*rho0

            m = bounce(l)

            dfeq = fequil_difference( &
                       l,rhocorner,ux,uy,0.0, &
                       inv1cs2,inv6cs6,ibgk)

            fedge = f(m,1,1,k) + dfeq

            f(l,0,0,k) = &
                 (1.0-win00)*fbase + &
                 win00*fedge

         endif

      endif


!=======================================================================
! Edge (i=0,j=ny+1)
!
! Adjacent completed face values:
!
!    i-face : f(l,0,ny,k)
!    j-face : f(l,1,ny+1,k)
!=======================================================================

      if (jN_is_phys) then

         fbase = ( ax*f(l,0,ny,k) + &
                   ay*f(l,1,ny+1,k) ) / aden

         f(l,0,ny+1,k) = fbase


!-----------------------------------------------------------------------
! Incoming through both faces:
!
!       cx > 0 and cy < 0
!-----------------------------------------------------------------------

         if (win0N > 0.0 .and. &
             cxs(l) > 0 .and. cys(l) < 0) then

            rholocal = 0.0

            do mm = 1,nl
               rholocal = rholocal + f(mm,1,ny,k)
            enddo

            rhocorner = rho_relax*rholocal + &
                        (1.0-rho_relax)*rho0

            m = bounce(l)

            dfeq = fequil_difference( &
                       l,rhocorner,ux,uy,0.0, &
                       inv1cs2,inv6cs6,ibgk)

            fedge = f(m,1,ny,k) + dfeq

            f(l,0,ny+1,k) = &
                 (1.0-win0N)*fbase + &
                 win0N*fedge

         endif

      endif


!=======================================================================
! Edge (i=nx+1,j=0)
!
! Adjacent completed face values:
!
!    i-face : f(l,nx+1,1,k)
!    j-face : f(l,nx,0,k)
!=======================================================================

      if (j0_is_phys) then

         fbase = ( ax*f(l,nx+1,1,k) + &
                   ay*f(l,nx,0,k) ) / aden

         f(l,nx+1,0,k) = fbase


!-----------------------------------------------------------------------
! Incoming through both faces:
!
!       cx < 0 and cy > 0
!-----------------------------------------------------------------------

         if (winN0 > 0.0 .and. &
             cxs(l) < 0 .and. cys(l) > 0) then

            rholocal = 0.0

            do mm = 1,nl
               rholocal = rholocal + f(mm,nx,1,k)
            enddo

            rhocorner = rho_relax*rholocal + &
                        (1.0-rho_relax)*rho0

            m = bounce(l)

            dfeq = fequil_difference( &
                       l,rhocorner,ux,uy,0.0, &
                       inv1cs2,inv6cs6,ibgk)

            fedge = f(m,nx,1,k) + dfeq

            f(l,nx+1,0,k) = &
                 (1.0-winN0)*fbase + &
                 winN0*fedge

         endif

      endif


!=======================================================================
! Edge (i=nx+1,j=ny+1)
!
! Adjacent completed face values:
!
!    i-face : f(l,nx+1,ny,k)
!    j-face : f(l,nx,ny+1,k)
!=======================================================================

      if (jN_is_phys) then

         fbase = ( ax*f(l,nx+1,ny,k) + &
                   ay*f(l,nx,ny+1,k) ) / aden

         f(l,nx+1,ny+1,k) = fbase


!-----------------------------------------------------------------------
! Incoming through both faces:
!
!       cx < 0 and cy < 0
!-----------------------------------------------------------------------

         if (winNN > 0.0 .and. &
             cxs(l) < 0 .and. cys(l) < 0) then

            rholocal = 0.0

            do mm = 1,nl
               rholocal = rholocal + f(mm,nx,ny,k)
            enddo

            rhocorner = rho_relax*rholocal + &
                        (1.0-rho_relax)*rho0

            m = bounce(l)

            dfeq = fequil_difference( &
                       l,rhocorner,ux,uy,0.0, &
                       inv1cs2,inv6cs6,ibgk)

            fedge = f(m,nx,ny,k) + dfeq

            f(l,nx+1,ny+1,k) = &
                 (1.0-winNN)*fbase + &
                 winNN*fedge

         endif

      endif


#ifndef _CUDA

   enddo
   enddo

!$OMP END PARALLEL DO

#endif


end subroutine boundary_inflow_edges_ij_kernel

end module m_boundary_inflow_edges_ij_kernel




!!!! module m_boundary_inflow_edges_ij_kernel
!!!! contains
!!!! 
!!!! #ifdef _CUDA
!!!! attributes(global) &
!!!! #endif
!!!! subroutine boundary_inflow_edges_ij_kernel( &
!!!!    f,uvel,udir,rho0,rho_relax,j0_is_phys,jN_is_phys)
!!!! 
!!!! !-----------------------------------------------------------------------
!!!! ! General inflow/zero-gradient treatment of the four true i-j corners
!!!! ! for k = 1..nz.
!!!! !
!!!! ! This routine is consistent with:
!!!! !
!!!! !    boundary_inflow_i_kernel
!!!! !    boundary_inflow_j_kernel
!!!! !    boundary_edges_mixed
!!!! !
!!!! ! and MUST be called after the i- and j-face boundary routines.
!!!! !
!!!! !
!!!! ! CORNER BASELINE
!!!! ! ---------------
!!!! !
!!!! ! Each corner is first constructed from the two already-completed
!!!! ! adjacent face ghost values:
!!!! !
!!!! !                  |ux| f_i-face + |uy| f_j-face
!!!! !       fbase = ----------------------------------
!!!! !                         |ux| + |uy|
!!!! !
!!!! ! This gives the correct limiting behaviour:
!!!! !
!!!! !    ux -> 0 : corner becomes continuation of the j-face condition
!!!! !    uy -> 0 : corner becomes continuation of the i-face condition
!!!! !
!!!! ! For example, at udir = 90 degrees:
!!!! !
!!!! !    f(:,0,0,k)       = f(:,1,0,k)
!!!! !    f(:,nx+1,0,k)    = f(:,nx,0,k)
!!!! !
!!!! ! so the lower corners are simply x-continuations of the j=0 inflow.
!!!! !
!!!! !
!!!! ! CORNER KRUGER RECONSTRUCTION
!!!! ! ----------------------------
!!!! !
!!!! ! If BOTH adjacent faces are inflow, populations which cross BOTH
!!!! ! boundaries are reconstructed using the Kruger condition.
!!!! !
!!!! ! The corresponding smooth corner weights are:
!!!! !
!!!! !    win00 = win_i0 * win_j0
!!!! !    win0N = win_i0 * win_jN
!!!! !    winN0 = win_iN * win_j0
!!!! !    winNN = win_iN * win_jN
!!!! !
!!!! ! Consequently the true two-face corner inflow smoothly disappears
!!!! ! whenever either normal velocity component approaches zero.
!!!! !
!!!! !
!!!! ! IMPORTANT:
!!!! !
!!!! ! Only populations having a nonzero incoming component through BOTH
!!!! ! boundaries are reconstructed at the corner:
!!!! !
!!!! !    (0,0)       : cx > 0, cy > 0
!!!! !    (0,ny+1)    : cx > 0, cy < 0
!!!! !    (nx+1,0)    : cx < 0, cy > 0
!!!! !    (nx+1,ny+1) : cx < 0, cy < 0
!!!! !
!!!! ! Populations with cx=0 or cy=0 are handled by the corresponding face
!!!! ! boundary and must NOT be treated as true corner-entering populations.
!!!! !
!!!! !
!!!! ! CUDA:
!!!! !
!!!! ! There are no thread-private arrays indexed by bounce(l).  The
!!!! ! required bounced population is evaluated directly, avoiding the
!!!! ! uninitialized-private-array problem in the previous implementation.
!!!! !
!!!! !
!!!! ! MPI:
!!!! !
!!!! ! j0_is_phys / jN_is_phys prevent modification of corners belonging
!!!! ! to internal MPI j-halo planes.
!!!! !
!!!! !-----------------------------------------------------------------------
!!!! 
!!!! #ifdef _CUDA
!!!!    use cudafor
!!!! #endif
!!!! 
!!!!    use mod_dimensions
!!!!    use mod_D3Q27setup, only : nl, weights, cxs, cys, bounce, cs2
!!!! 
!!!!    implicit none
!!!! 
!!!!    real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
!!!!    real, intent(in)    :: uvel(nz)
!!!! 
!!!!    real, value    :: udir
!!!!    real, value    :: rho0
!!!!    real, value    :: rho_relax
!!!! 
!!!!    logical, value :: j0_is_phys
!!!!    logical, value :: jN_is_phys
!!!! 
!!!! #ifdef _CUDA
!!!!    attributes(device) :: f
!!!!    attributes(device) :: uvel
!!!! #endif
!!!! 
!!!!    integer :: k,l,m,mm
!!!! 
!!!!    real, parameter :: pi          = acos(-1.0)
!!!!    real, parameter :: blend_width = 0.05
!!!! 
!!!!    real :: uxdir,uydir
!!!!    real :: ax,ay,aden
!!!! 
!!!!    real :: uu
!!!!    real :: s
!!!! 
!!!!    real :: win_i0,win_iN
!!!!    real :: win_j0,win_jN
!!!! 
!!!!    real :: win00,win0N
!!!!    real :: winN0,winNN
!!!! 
!!!!    real :: fbase
!!!!    real :: fkruger
!!!! 
!!!!    real :: rholocal
!!!!    real :: rhocorner
!!!! 
!!!!    real :: momentum_correction
!!!! 
!!!! 
!!!! !=======================================================================
!!!! ! CUDA / CPU indexing
!!!! !=======================================================================
!!!! 
!!!! #ifdef _CUDA
!!!! 
!!!!    l = threadIdx%x + (blockIdx%x-1)*blockDim%x
!!!!    k = threadIdx%y + (blockIdx%y-1)*blockDim%y
!!!! 
!!!!    if (l < 1 .or. l > nl) return
!!!!    if (k < 1 .or. k > nz) return
!!!! 
!!!! #else
!!!! 
!!!! !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                          &
!!!! !$OMP& PRIVATE(l,k,m,mm,uxdir,uydir,ax,ay,aden,uu,s,                &
!!!! !$OMP&         win_i0,win_iN,win_j0,win_jN,                         &
!!!! !$OMP&         win00,win0N,winN0,winNN,                             &
!!!! !$OMP&         fbase,fkruger,rholocal,rhocorner,                    &
!!!! !$OMP&         momentum_correction)                                  &
!!!! !$OMP& SHARED(f,uvel,udir,rho0,rho_relax,j0_is_phys,jN_is_phys)
!!!! 
!!!!    do k = 1,nz
!!!!    do l = 1,nl
!!!! 
!!!! #endif
!!!! 
!!!! 
!!!! !=======================================================================
!!!! ! Prescribed velocity direction
!!!! !=======================================================================
!!!! 
!!!!       uxdir = cos(udir*pi/180.0)
!!!!       uydir = sin(udir*pi/180.0)
!!!! 
!!!!       uu = uvel(k)
!!!! 
!!!!       ax = abs(uxdir)
!!!!       ay = abs(uydir)
!!!! 
!!!!       ! uxdir and uydir form a unit vector, so aden is always > 0.
!!!!       aden = ax + ay
!!!! 
!!!! 
!!!! !=======================================================================
!!!! ! Smooth inflow weights.
!!!! !
!!!! ! These must use exactly the same definition and blend_width as the
!!!! ! corresponding i- and j-face kernels.
!!!! !=======================================================================
!!!! 
!!!!       !---------------------------------------------------------------
!!!!       ! i = 0 inflow
!!!!       !---------------------------------------------------------------
!!!! 
!!!!       s = min(1.0,max(0.0,uxdir/blend_width))
!!!! 
!!!!       win_i0 = s*s*(3.0-2.0*s)
!!!! 
!!!! 
!!!!       !---------------------------------------------------------------
!!!!       ! i = nx+1 inflow
!!!!       !---------------------------------------------------------------
!!!! 
!!!!       s = min(1.0,max(0.0,-uxdir/blend_width))
!!!! 
!!!!       win_iN = s*s*(3.0-2.0*s)
!!!! 
!!!! 
!!!!       !---------------------------------------------------------------
!!!!       ! j = 0 inflow
!!!!       !---------------------------------------------------------------
!!!! 
!!!!       s = min(1.0,max(0.0,uydir/blend_width))
!!!! 
!!!!       win_j0 = s*s*(3.0-2.0*s)
!!!! 
!!!! 
!!!!       !---------------------------------------------------------------
!!!!       ! j = ny+1 inflow
!!!!       !---------------------------------------------------------------
!!!! 
!!!!       s = min(1.0,max(0.0,-uydir/blend_width))
!!!! 
!!!!       win_jN = s*s*(3.0-2.0*s)
!!!! 
!!!! 
!!!! !=======================================================================
!!!! ! True two-face corner inflow weights.
!!!! !=======================================================================
!!!! 
!!!!       win00 = win_i0*win_j0
!!!!       win0N = win_i0*win_jN
!!!! 
!!!!       winN0 = win_iN*win_j0
!!!!       winNN = win_iN*win_jN
!!!! 
!!!! 
!!!! !=======================================================================
!!!! ! Corner (i=0,j=0)
!!!! !
!!!! ! Adjacent completed face values:
!!!! !
!!!! !    i-face : f(l,0,1,k)
!!!! !    j-face : f(l,1,0,k)
!!!! !
!!!! ! At ux=0 this becomes exactly f(l,1,0,k).
!!!! ! At uy=0 this becomes exactly f(l,0,1,k).
!!!! !=======================================================================
!!!! 
!!!!       if (j0_is_phys) then
!!!! 
!!!!          fbase = ( ax*f(l,0,1,k) + &
!!!!                    ay*f(l,1,0,k) ) / aden
!!!! 
!!!!          f(l,0,0,k) = fbase
!!!! 
!!!! 
!!!!          !------------------------------------------------------------
!!!!          ! True corner-entering populations:
!!!!          !
!!!!          !       cx > 0 and cy > 0
!!!!          !------------------------------------------------------------
!!!! 
!!!!          if (win00 > 0.0) then
!!!! 
!!!!             if (cxs(l) > 0 .and. cys(l) > 0) then
!!!! 
!!!!                rholocal = 0.0
!!!! 
!!!!                do mm = 1,nl
!!!!                   rholocal = rholocal + f(mm,1,1,k)
!!!!                enddo
!!!! 
!!!!                rhocorner = rho_relax*rholocal + &
!!!!                            (1.0-rho_relax)*rho0
!!!! 
!!!! 
!!!!                !-----------------------------------------------------
!!!!                ! Direct evaluation of the bounced Kruger population.
!!!!                !
!!!!                ! Since m=bounce(l):
!!!!                !
!!!!                !    c(m) = -c(l)
!!!!                !
!!!!                ! giving
!!!!                !
!!!!                ! fkruger =
!!!!                !    f(m) + 2*w(l)*rho*c(l).u/cs2
!!!!                !
!!!!                !-----------------------------------------------------
!!!! 
!!!!                m = bounce(l)
!!!! 
!!!!                momentum_correction =                         &
!!!!                     2.0*weights(l)*rhocorner*uu             &
!!!!                     *( real(cxs(l))*uxdir +                 &
!!!!                        real(cys(l))*uydir ) / cs2
!!!! 
!!!!                fkruger = f(m,1,1,k) + momentum_correction
!!!! 
!!!! 
!!!!                f(l,0,0,k) = &
!!!!                     (1.0-win00)*fbase + &
!!!!                     win00*fkruger
!!!! 
!!!!             endif
!!!! 
!!!!          endif
!!!! 
!!!!       endif
!!!! 
!!!! 
!!!! !=======================================================================
!!!! ! Corner (i=0,j=ny+1)
!!!! !
!!!! ! Adjacent completed face values:
!!!! !
!!!! !    i-face : f(l,0,ny,k)
!!!! !    j-face : f(l,1,ny+1,k)
!!!! !=======================================================================
!!!! 
!!!!       if (jN_is_phys) then
!!!! 
!!!!          fbase = ( ax*f(l,0,ny,k) + &
!!!!                    ay*f(l,1,ny+1,k) ) / aden
!!!! 
!!!!          f(l,0,ny+1,k) = fbase
!!!! 
!!!! 
!!!!          !------------------------------------------------------------
!!!!          ! True corner-entering populations:
!!!!          !
!!!!          !       cx > 0 and cy < 0
!!!!          !------------------------------------------------------------
!!!! 
!!!!          if (win0N > 0.0) then
!!!! 
!!!!             if (cxs(l) > 0 .and. cys(l) < 0) then
!!!! 
!!!!                rholocal = 0.0
!!!! 
!!!!                do mm = 1,nl
!!!!                   rholocal = rholocal + f(mm,1,ny,k)
!!!!                enddo
!!!! 
!!!!                rhocorner = rho_relax*rholocal + &
!!!!                            (1.0-rho_relax)*rho0
!!!! 
!!!!                m = bounce(l)
!!!! 
!!!!                momentum_correction =                         &
!!!!                     2.0*weights(l)*rhocorner*uu             &
!!!!                     *( real(cxs(l))*uxdir +                 &
!!!!                        real(cys(l))*uydir ) / cs2
!!!! 
!!!!                fkruger = f(m,1,ny,k) + momentum_correction
!!!! 
!!!! 
!!!!                f(l,0,ny+1,k) = &
!!!!                     (1.0-win0N)*fbase + &
!!!!                     win0N*fkruger
!!!! 
!!!!             endif
!!!! 
!!!!          endif
!!!! 
!!!!       endif
!!!! 
!!!! 
!!!! !=======================================================================
!!!! ! Corner (i=nx+1,j=0)
!!!! !
!!!! ! Adjacent completed face values:
!!!! !
!!!! !    i-face : f(l,nx+1,1,k)
!!!! !    j-face : f(l,nx,0,k)
!!!! !=======================================================================
!!!! 
!!!!       if (j0_is_phys) then
!!!! 
!!!!          fbase = ( ax*f(l,nx+1,1,k) + &
!!!!                    ay*f(l,nx,0,k) ) / aden
!!!! 
!!!!          f(l,nx+1,0,k) = fbase
!!!! 
!!!! 
!!!!          !------------------------------------------------------------
!!!!          ! True corner-entering populations:
!!!!          !
!!!!          !       cx < 0 and cy > 0
!!!!          !------------------------------------------------------------
!!!! 
!!!!          if (winN0 > 0.0) then
!!!! 
!!!!             if (cxs(l) < 0 .and. cys(l) > 0) then
!!!! 
!!!!                rholocal = 0.0
!!!! 
!!!!                do mm = 1,nl
!!!!                   rholocal = rholocal + f(mm,nx,1,k)
!!!!                enddo
!!!! 
!!!!                rhocorner = rho_relax*rholocal + &
!!!!                            (1.0-rho_relax)*rho0
!!!! 
!!!!                m = bounce(l)
!!!! 
!!!!                momentum_correction =                         &
!!!!                     2.0*weights(l)*rhocorner*uu             &
!!!!                     *( real(cxs(l))*uxdir +                 &
!!!!                        real(cys(l))*uydir ) / cs2
!!!! 
!!!!                fkruger = f(m,nx,1,k) + momentum_correction
!!!! 
!!!! 
!!!!                f(l,nx+1,0,k) = &
!!!!                     (1.0-winN0)*fbase + &
!!!!                     winN0*fkruger
!!!! 
!!!!             endif
!!!! 
!!!!          endif
!!!! 
!!!!       endif
!!!! 
!!!! 
!!!! !=======================================================================
!!!! ! Corner (i=nx+1,j=ny+1)
!!!! !
!!!! ! Adjacent completed face values:
!!!! !
!!!! !    i-face : f(l,nx+1,ny,k)
!!!! !    j-face : f(l,nx,ny+1,k)
!!!! !=======================================================================
!!!! 
!!!!       if (jN_is_phys) then
!!!! 
!!!!          fbase = ( ax*f(l,nx+1,ny,k) + &
!!!!                    ay*f(l,nx,ny+1,k) ) / aden
!!!! 
!!!!          f(l,nx+1,ny+1,k) = fbase
!!!! 
!!!! 
!!!!          !------------------------------------------------------------
!!!!          ! True corner-entering populations:
!!!!          !
!!!!          !       cx < 0 and cy < 0
!!!!          !------------------------------------------------------------
!!!! 
!!!!          if (winNN > 0.0) then
!!!! 
!!!!             if (cxs(l) < 0 .and. cys(l) < 0) then
!!!! 
!!!!                rholocal = 0.0
!!!! 
!!!!                do mm = 1,nl
!!!!                   rholocal = rholocal + f(mm,nx,ny,k)
!!!!                enddo
!!!! 
!!!!                rhocorner = rho_relax*rholocal + &
!!!!                            (1.0-rho_relax)*rho0
!!!! 
!!!!                m = bounce(l)
!!!! 
!!!!                momentum_correction =                         &
!!!!                     2.0*weights(l)*rhocorner*uu             &
!!!!                     *( real(cxs(l))*uxdir +                 &
!!!!                        real(cys(l))*uydir ) / cs2
!!!! 
!!!!                fkruger = f(m,nx,ny,k) + momentum_correction
!!!! 
!!!! 
!!!!                f(l,nx+1,ny+1,k) = &
!!!!                     (1.0-winNN)*fbase + &
!!!!                     winNN*fkruger
!!!! 
!!!!             endif
!!!! 
!!!!          endif
!!!! 
!!!!       endif
!!!! 
!!!! 
!!!! #ifndef _CUDA
!!!! 
!!!!    enddo
!!!!    enddo
!!!! 
!!!! !$OMP END PARALLEL DO
!!!! 
!!!! #endif
!!!! 
!!!! 
!!!! end subroutine boundary_inflow_edges_ij_kernel
!!!! 
!!!! end module m_boundary_inflow_edges_ij_kernel
!!!!New Claude
!!!module m_boundary_inflow_edges_ij_kernel
!!!contains
!!!
!!!#ifdef _CUDA
!!!   attributes(global) &
!!!#endif
!!!subroutine boundary_inflow_edges_ij_kernel(f,uvel,udir,rho0,rho_relax,j0_is_phys,jN_is_phys)
!!!
!!!!-----------------------------------------------------------------------
!!!! General inflow/zero-gradient boundary condition at the four true
!!!! i-j corners (k = 1..nz), consistent with boundary_inflow_i_kernel
!!!! and boundary_inflow_j_kernel.
!!!!
!!!! Each corner's Kruger weight is the PRODUCT of the corresponding
!!!! i-face and j-face smoothstep weights:
!!!!
!!!!     win_00 = win_i0 * win_j0   (corner i=0,    j=0   )
!!!!     win_0N = win_i0 * win_jN   (corner i=0,    j=ny+1)
!!!!     win_N0 = win_iN * win_j0   (corner i=nx+1, j=0   )
!!!!     win_NN = win_iN * win_jN   (corner i=nx+1, j=ny+1)
!!!!
!!!! so a given corner is only fully "Kruger inflow" when BOTH adjacent
!!!! faces are simultaneously inflow-aligned beyond blend_width; if
!!!! either uxdir or uydir crosses toward tangential/outflow for that
!!!! corner, the weight smoothly collapses to 0 and the corner falls
!!!! back to zero-normal-gradient (direct copy from the nearest interior
!!!! corner-adjacent node) - no history-dependent convective term, same
!!!! as the face kernels.
!!!!
!!!! j0_is_phys / jN_is_phys gate the two corners on each j-plane the
!!!! same way as boundary_inflow_j_kernel, for MPI interfaces. i is never
!!!! MPI-decomposed here, so both i-side corners are always physical.
!!!!
!!!! blend_width matches the face kernels' convention: expressed in
!!!! direction-cosine units, 0.05 corresponds to about +/- 2.87 degrees.
!!!!-----------------------------------------------------------------------
!!!
!!!#ifdef _CUDA
!!!   use cudafor
!!!#endif
!!!
!!!   use mod_dimensions
!!!   use mod_D3Q27setup, only : nl, weights, cxs, cys, bounce, cs2
!!!
!!!   implicit none
!!!
!!!   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
!!!   real, intent(in)    :: uvel(nz)
!!!
!!!   real, value    :: udir
!!!   real, value    :: rho0
!!!   real, value    :: rho_relax
!!!   logical, value :: j0_is_phys
!!!   logical, value :: jN_is_phys
!!!
!!!   integer :: k,l
!!!
!!!   real, parameter :: pi          = acos(-1.0)
!!!   real, parameter :: blend_width = 0.05
!!!
!!!   real :: wl, cxl, cyl
!!!   real :: uu
!!!   real :: uxdir, uydir
!!!   real :: s
!!!   real :: win_i0, win_iN, win_j0, win_jN
!!!   real :: win00, win0N, winN0, winNN
!!!
!!!   real :: rholocal, rhocorner
!!!   real :: fraw00(nl), fraw0N(nl), frawN0(nl), frawNN(nl)
!!!   real :: fkr00(nl),  fkr0N(nl),  fkrN0(nl),  fkrNN(nl)
!!!
!!!#ifdef _CUDA
!!!   l = threadIdx%x + (blockIdx%x-1)*blockDim%x
!!!   k = threadIdx%y + (blockIdx%y-1)*blockDim%y
!!!   if (l < 1 .or. l > nl) return
!!!   if (k < 1 .or. k > nz) return
!!!#else
!!!!$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                                    &
!!!!$OMP& PRIVATE(l,k,wl,cxl,cyl,uu,uxdir,uydir,s,                                &
!!!!$OMP&         win_i0,win_iN,win_j0,win_jN,win00,win0N,winN0,winNN,            &
!!!!$OMP&         rholocal,rhocorner,                                             &
!!!!$OMP&         fraw00,fraw0N,frawN0,frawNN,fkr00,fkr0N,fkrN0,fkrNN)            &
!!!!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,j0_is_phys,jN_is_phys)
!!!   do k = 1,nz
!!!   do l = 1,nl
!!!#endif
!!!
!!!      uxdir = cos(udir*pi/180.0)
!!!      uydir = sin(udir*pi/180.0)
!!!      uu    = uvel(k)
!!!
!!!      !--------------------------------------------------------------
!!!      ! Per-face smoothstep weights, same convention as the face
!!!      ! kernels' win0/winN.
!!!      !--------------------------------------------------------------
!!!      s = min(1.0,max(0.0, uxdir/blend_width));  win_i0 = s*s*(3.0-2.0*s)
!!!      s = min(1.0,max(0.0,-uxdir/blend_width));  win_iN = s*s*(3.0-2.0*s)
!!!      s = min(1.0,max(0.0, uydir/blend_width));  win_j0 = s*s*(3.0-2.0*s)
!!!      s = min(1.0,max(0.0,-uydir/blend_width));  win_jN = s*s*(3.0-2.0*s)
!!!
!!!      win00 = win_i0*win_j0   ! corner (0,0)
!!!      win0N = win_i0*win_jN   ! corner (0,ny+1)
!!!      winN0 = win_iN*win_j0   ! corner (nx+1,0)
!!!      winNN = win_iN*win_jN   ! corner (nx+1,ny+1)
!!!
!!!      !--------------------------------------------------------------
!!!      ! Local density at each corner's nearest interior node.
!!!      !--------------------------------------------------------------
!!!      rholocal = 0.0
!!!      do wl = 1,1  ! placeholder loop removed below; see explicit sums
!!!      enddo
!!!
!!!      ! Corner (0,0), reference node (1,1,k)
!!!      rholocal = 0.0
!!!      block
!!!         integer :: mm
!!!         do mm = 1,nl
!!!            rholocal = rholocal + f(mm,1,1,k)
!!!         enddo
!!!      end block
!!!      rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
!!!      wl  = weights(l); cxl = real(cxs(l)); cyl = real(cys(l))
!!!      fraw00(l) = f(l,1,1,k) - 2.0*wl*rhocorner*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
!!!
!!!      ! Corner (0,ny+1), reference node (1,ny,k)
!!!      rholocal = 0.0
!!!      block
!!!         integer :: mm
!!!         do mm = 1,nl
!!!            rholocal = rholocal + f(mm,1,ny,k)
!!!         enddo
!!!      end block
!!!      rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
!!!      fraw0N(l) = f(l,1,ny,k) - 2.0*wl*rhocorner*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
!!!
!!!      ! Corner (nx+1,0), reference node (nx,1,k)
!!!      rholocal = 0.0
!!!      block
!!!         integer :: mm
!!!         do mm = 1,nl
!!!            rholocal = rholocal + f(mm,nx,1,k)
!!!         enddo
!!!      end block
!!!      rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
!!!      frawN0(l) = f(l,nx,1,k) - 2.0*wl*rhocorner*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
!!!
!!!      ! Corner (nx+1,ny+1), reference node (nx,ny,k)
!!!      rholocal = 0.0
!!!      block
!!!         integer :: mm
!!!         do mm = 1,nl
!!!            rholocal = rholocal + f(mm,nx,ny,k)
!!!         enddo
!!!      end block
!!!      rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
!!!      frawNN(l) = f(l,nx,ny,k) - 2.0*wl*rhocorner*(cxl*uu*uxdir + cyl*uu*uydir)/cs2
!!!
!!!      !--------------------------------------------------------------
!!!      ! Bounce mapping: a population "enters" a given corner when its
!!!      ! direction has the correct sign in BOTH i and j simultaneously
!!!      ! (diagonal or single-axis-into-the-corner directions alike).
!!!      !--------------------------------------------------------------
!!!      if (cxs(l) >= 0 .and. cys(l) >= 0) then
!!!         fkr00(l) = fraw00(bounce(l))
!!!      else
!!!         fkr00(l) = fraw00(l)
!!!      endif
!!!
!!!      if (cxs(l) >= 0 .and. cys(l) <= 0) then
!!!         fkr0N(l) = fraw0N(bounce(l))
!!!      else
!!!         fkr0N(l) = fraw0N(l)
!!!      endif
!!!
!!!      if (cxs(l) <= 0 .and. cys(l) >= 0) then
!!!         fkrN0(l) = frawN0(bounce(l))
!!!      else
!!!         fkrN0(l) = frawN0(l)
!!!      endif
!!!
!!!      if (cxs(l) <= 0 .and. cys(l) <= 0) then
!!!         fkrNN(l) = frawNN(bounce(l))
!!!      else
!!!         fkrNN(l) = frawNN(l)
!!!      endif
!!!
!!!      !--------------------------------------------------------------
!!!      ! Final blend: win=1 -> Kruger corner inflow; win=0 -> plain
!!!      ! zero-normal-gradient copy from the nearest interior node.
!!!      ! Gated by j0_is_phys/jN_is_phys for MPI interfaces.
!!!      !--------------------------------------------------------------
!!!      if (j0_is_phys) then
!!!         f(l,0,0,k)    = win00*fkr00(l) + (1.0-win00)*f(l,1,1,k)
!!!         f(l,nx+1,0,k) = winN0*fkrN0(l) + (1.0-winN0)*f(l,nx,1,k)
!!!      endif
!!!
!!!      if (jN_is_phys) then
!!!         f(l,0,ny+1,k)    = win0N*fkr0N(l) + (1.0-win0N)*f(l,1,ny,k)
!!!         f(l,nx+1,ny+1,k) = winNN*fkrNN(l) + (1.0-winNN)*f(l,nx,ny,k)
!!!      endif
!!!
!!!#ifndef _CUDA
!!!   enddo
!!!   enddo
!!!!$OMP END PARALLEL DO
!!!#endif
!!!
!!!end subroutine boundary_inflow_edges_ij_kernel
!!!
!!!end module m_boundary_inflow_edges_ij_kernel
!!!



! New ChatGPT
!! module m_boundary_inflow_edges_ij_kernel
!! contains
!! 
!! #ifdef _CUDA
!! attributes(global) &
!! #endif
!! subroutine boundary_inflow_edges_ij_kernel(f,uvel,udir,rho0,rho_relax)
!! 
!! !-----------------------------------------------------------------------
!! ! Treatment of the four i-j corner/edge ghost lines.
!! !
!! ! This routine must be called AFTER the i- and j-face boundary routines,
!! ! because it uses the already updated face ghost values.
!! !
!! ! The corner treatment is consistent with the smooth face treatment:
!! !
!! !   1. A baseline corner value is formed from the two adjacent face
!! !      ghost values.
!! !
!! !   2. The relative weights of the two faces are proportional to the
!! !      magnitude of the corresponding normal velocity components:
!! !
!! !          |ux| and |uy|.
!! !
!! !      Therefore:
!! !
!! !        ux -> 0 : corner becomes an extrapolation of the j-face BC
!! !        uy -> 0 : corner becomes an extrapolation of the i-face BC
!! !
!! !   3. If BOTH adjacent faces are inflow, the populations entering
!! !      through both boundaries are reconstructed with the Kruger
!! !      condition.
!! !
!! !   4. This corner reconstruction is switched on smoothly using the
!! !      product of the two corresponding face-inflow weights.
!! !
!! ! Thus, for example, around udir = 90 deg:
!! !
!! !      ux -> 0
!! !      uy > 0
!! !
!! ! the two lower corners smoothly approach
!! !
!! !      f(:,0,0)       -> f(:,1,0)
!! !      f(:,nx+1,0)    -> f(:,nx,0)
!! !
!! ! i.e. simple x-extrapolation of the j=0 inflow boundary.
!! !
!! ! No old corner ghost value enters the calculation, eliminating the
!! ! history-dependent convective corner condition.
!! !
!! !-----------------------------------------------------------------------
!! 
!! #ifdef _CUDA
!!    use cudafor
!! #endif
!! 
!!    use mod_dimensions
!!    use mod_D3Q27setup, only : nl,weights,cxs,cys,bounce,cs2
!! 
!!    implicit none
!! 
!!    real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
!!    real, intent(in)    :: uvel(nz)
!! 
!!    real, value :: udir
!!    real, value :: rho0
!!    real, value :: rho_relax
!! 
!! #ifdef _CUDA
!!    attributes(device) :: f
!!    attributes(device) :: uvel
!! #endif
!! 
!!    integer :: k,l,m,mm
!! 
!!    real, parameter :: pi = acos(-1.0)
!! 
!!    ! Must normally be the same value as in the i- and j-face kernels.
!!    real, parameter :: blend_width = 0.05
!! 
!!    real :: uxdir,uydir
!!    real :: ax,ay,aden
!! 
!!    real :: uu
!! 
!!    real :: sx0,sxN
!!    real :: sy0,syN
!! 
!!    real :: wx0,wxN
!!    real :: wy0,wyN
!! 
!!    real :: wc00,wc0N,wcN0,wcNN
!! 
!!    real :: fbase
!!    real :: fkruger
!! 
!!    real :: rholocal
!!    real :: rhocorner
!! 
!!    real :: momentum_correction
!! 
!! 
!! !=======================================================================
!! ! CUDA / CPU indexing
!! !=======================================================================
!! 
!! #ifdef _CUDA
!! 
!!    l = threadIdx%x + (blockIdx%x-1)*blockDim%x
!!    k = threadIdx%y + (blockIdx%y-1)*blockDim%y
!! 
!!    if (l < 1 .or. l > nl) return
!!    if (k < 1 .or. k > nz) return
!! 
!! #else
!! 
!! !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                           &
!! !$OMP PRIVATE(l,k,m,mm,uxdir,uydir,ax,ay,aden,uu,                    &
!! !$OMP         sx0,sxN,sy0,syN,wx0,wxN,wy0,wyN,                      &
!! !$OMP         wc00,wc0N,wcN0,wcNN,                                  &
!! !$OMP         fbase,fkruger,rholocal,rhocorner,                      &
!! !$OMP         momentum_correction)                                   &
!! !$OMP SHARED(f,uvel,udir,rho0,rho_relax)
!! 
!!    do k = 1,nz
!!    do l = 1,nl
!! 
!! #endif
!! 
!! 
!! !=======================================================================
!! ! Wind direction
!! !=======================================================================
!! 
!!       uxdir = cos(udir*pi/180.0)
!!       uydir = sin(udir*pi/180.0)
!! 
!!       uu = uvel(k)
!! 
!!       ax = abs(uxdir)
!!       ay = abs(uydir)
!! 
!!       ! Since uxdir^2 + uydir^2 = 1 this denominator cannot vanish.
!!       aden = ax + ay
!! 
!! 
!! !=======================================================================
!! ! Smooth face-inflow weights.
!! !
!! ! These are exactly analogous to those used in the two face kernels.
!! !
!! ! wx0 : inflow through i=0
!! ! wxN : inflow through i=nx+1
!! !
!! ! wy0 : inflow through j=0
!! ! wyN : inflow through j=ny+1
!! !=======================================================================
!! 
!!       sx0 = min(1.0,max(0.0, uxdir/blend_width))
!!       sxN = min(1.0,max(0.0,-uxdir/blend_width))
!! 
!!       sy0 = min(1.0,max(0.0, uydir/blend_width))
!!       syN = min(1.0,max(0.0,-uydir/blend_width))
!! 
!! 
!!       wx0 = sx0*sx0*(3.0-2.0*sx0)
!!       wxN = sxN*sxN*(3.0-2.0*sxN)
!! 
!!       wy0 = sy0*sy0*(3.0-2.0*sy0)
!!       wyN = syN*syN*(3.0-2.0*syN)
!! 
!! 
!! !=======================================================================
!! ! Corner inflow weights.
!! !
!! ! A genuine corner inflow exists only when BOTH adjacent faces are
!! ! inflow.
!! !=======================================================================
!! 
!!       wc00 = wx0*wy0
!!       wc0N = wx0*wyN
!!       wcN0 = wxN*wy0
!!       wcNN = wxN*wyN
!! 
!! 
!! !=======================================================================
!! ! Corner (0,0)
!! !
!! ! Baseline:
!! !
!! !    x-face contribution : f(l,0,1,k)
!! !    j-face contribution : f(l,1,0,k)
!! !
!! ! At udir=90:
!! !
!! !    ax=0, ay=1
!! !
!! ! so this becomes exactly
!! !
!! !    f(l,1,0,k)
!! !
!! ! which is the already constructed j=0 inflow ghost value.
!! !=======================================================================
!! 
!!       fbase = ( ax*f(l,0,1,k) + &
!!                 ay*f(l,1,0,k) ) / aden
!! 
!!       f(l,0,0,k) = fbase
!! 
!! 
!!       !----------------------------------------------------------------
!!       ! If both i=0 and j=0 are inflow, reconstruct populations which
!!       ! enter the domain through BOTH faces:
!!       !
!!       !          cx > 0, cy > 0
!!       !----------------------------------------------------------------
!! 
!!       if (wc00 > 0.0) then
!! 
!!          if (cxs(l) > 0 .and. cys(l) > 0) then
!! 
!!             rholocal = 0.0
!! 
!!             do mm = 1,nl
!!                rholocal = rholocal + f(mm,1,1,k)
!!             enddo
!! 
!!             rhocorner = rho_relax*rholocal + &
!!                         (1.0-rho_relax)*rho0
!! 
!!             m = bounce(l)
!! 
!!             momentum_correction =                         &
!!                  2.0*weights(l)*rhocorner*uu             &
!!                  *(real(cxs(l))*uxdir +                  &
!!                    real(cys(l))*uydir)/cs2
!! 
!!             fkruger = f(m,1,1,k) + momentum_correction
!! 
!!             f(l,0,0,k) = (1.0-wc00)*fbase + &
!!                           wc00*fkruger
!! 
!!          endif
!! 
!!       endif
!! 
!! 
!! !=======================================================================
!! ! Corner (0,ny+1)
!! !=======================================================================
!! 
!!       fbase = ( ax*f(l,0,ny,k) + &
!!                 ay*f(l,1,ny+1,k) ) / aden
!! 
!!       f(l,0,ny+1,k) = fbase
!! 
!! 
!!       !----------------------------------------------------------------
!!       ! Inflow through both faces:
!!       !
!!       !        i=0       -> cx > 0
!!       !        j=ny+1    -> cy < 0
!!       !----------------------------------------------------------------
!! 
!!       if (wc0N > 0.0) then
!! 
!!          if (cxs(l) > 0 .and. cys(l) < 0) then
!! 
!!             rholocal = 0.0
!! 
!!             do mm = 1,nl
!!                rholocal = rholocal + f(mm,1,ny,k)
!!             enddo
!! 
!!             rhocorner = rho_relax*rholocal + &
!!                         (1.0-rho_relax)*rho0
!! 
!!             m = bounce(l)
!! 
!!             momentum_correction =                         &
!!                  2.0*weights(l)*rhocorner*uu             &
!!                  *(real(cxs(l))*uxdir +                  &
!!                    real(cys(l))*uydir)/cs2
!! 
!!             fkruger = f(m,1,ny,k) + momentum_correction
!! 
!!             f(l,0,ny+1,k) = (1.0-wc0N)*fbase + &
!!                              wc0N*fkruger
!! 
!!          endif
!! 
!!       endif
!! 
!! 
!! !=======================================================================
!! ! Corner (nx+1,0)
!! !=======================================================================
!! 
!!       fbase = ( ax*f(l,nx+1,1,k) + &
!!                 ay*f(l,nx,0,k) ) / aden
!! 
!!       f(l,nx+1,0,k) = fbase
!! 
!! 
!!       !----------------------------------------------------------------
!!       ! Inflow through both faces:
!!       !
!!       !        i=nx+1    -> cx < 0
!!       !        j=0       -> cy > 0
!!       !----------------------------------------------------------------
!! 
!!       if (wcN0 > 0.0) then
!! 
!!          if (cxs(l) < 0 .and. cys(l) > 0) then
!! 
!!             rholocal = 0.0
!! 
!!             do mm = 1,nl
!!                rholocal = rholocal + f(mm,nx,1,k)
!!             enddo
!! 
!!             rhocorner = rho_relax*rholocal + &
!!                         (1.0-rho_relax)*rho0
!! 
!!             m = bounce(l)
!! 
!!             momentum_correction =                         &
!!                  2.0*weights(l)*rhocorner*uu             &
!!                  *(real(cxs(l))*uxdir +                  &
!!                    real(cys(l))*uydir)/cs2
!! 
!!             fkruger = f(m,nx,1,k) + momentum_correction
!! 
!!             f(l,nx+1,0,k) = (1.0-wcN0)*fbase + &
!!                              wcN0*fkruger
!! 
!!          endif
!! 
!!       endif
!! 
!! 
!! !=======================================================================
!! ! Corner (nx+1,ny+1)
!! !=======================================================================
!! 
!!       fbase = ( ax*f(l,nx+1,ny,k) + &
!!                 ay*f(l,nx,ny+1,k) ) / aden
!! 
!!       f(l,nx+1,ny+1,k) = fbase
!! 
!! 
!!       !----------------------------------------------------------------
!!       ! Inflow through both faces:
!!       !
!!       !        i=nx+1    -> cx < 0
!!       !        j=ny+1    -> cy < 0
!!       !----------------------------------------------------------------
!! 
!!       if (wcNN > 0.0) then
!! 
!!          if (cxs(l) < 0 .and. cys(l) < 0) then
!! 
!!             rholocal = 0.0
!! 
!!             do mm = 1,nl
!!                rholocal = rholocal + f(mm,nx,ny,k)
!!             enddo
!! 
!!             rhocorner = rho_relax*rholocal + &
!!                         (1.0-rho_relax)*rho0
!! 
!!             m = bounce(l)
!! 
!!             momentum_correction =                         &
!!                  2.0*weights(l)*rhocorner*uu             &
!!                  *(real(cxs(l))*uxdir +                  &
!!                    real(cys(l))*uydir)/cs2
!! 
!!             fkruger = f(m,nx,ny,k) + momentum_correction
!! 
!!             f(l,nx+1,ny+1,k) = (1.0-wcNN)*fbase + &
!!                                 wcNN*fkruger
!! 
!!          endif
!! 
!!       endif
!! 
!! 
!! #ifndef _CUDA
!! 
!!    enddo
!!    enddo
!! 
!! !$OMP END PARALLEL DO
!! 
!! #endif
!! 
!! 
!! end subroutine boundary_inflow_edges_ij_kernel
!! 
!! end module m_boundary_inflow_edges_ij_kernel

!old Claude
!! module m_boundary_inflow_edges_ij_kernel
!! contains
!! #ifdef _CUDA
!! attributes(global) &
!! #endif
!! subroutine boundary_inflow_edges_ij_kernel(f,uvel,udir,rho0,rho_relax)
!! 
!! #ifdef _CUDA
!!    use cudafor
!! #endif
!!    use mod_dimensions
!!    use mod_D3Q27setup, only : nl, weights, cxs, cys, bounce, cs2
!! 
!!    implicit none
!! 
!!    real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
!!    real, intent(in)    :: uvel(nz)
!! 
!!    real, value :: udir
!!    real, value :: rho0
!!    real, value :: rho_relax
!!    real        :: uconv
!! 
!!    real, parameter :: udir_tol = 10.0*epsilon(1.0)
!!    real, parameter :: pi       = acos(-1.0)
!!    real, parameter :: uconv_min_frac = 0.2   ! same floor fraction as face kernels
!! 
!! #ifdef _CUDA
!!    attributes(device) :: f
!!    attributes(device) :: uvel
!! #endif
!! 
!!    integer :: k,l,m,mm
!!    real :: uxdir,uydir
!!    real :: uu,rhocorner,rholocal
!!    real :: momentum_correction
!!    real :: cconv, invden
!!    real :: uvel_ref   ! reference bulk speed for the floor (see note below)
!! 
!! #ifdef _CUDA
!!    l = threadIdx%x + (blockIdx%x-1)*blockDim%x
!!    k = threadIdx%y + (blockIdx%y-1)*blockDim%y
!!    if (l < 1 .or. l > nl) return
!!    if (k < 1 .or. k > nz) return
!! #else
!! !$OMP PARALLEL DO COLLAPSE(2) DEFAULT(NONE)                                        &
!! !$OMP PRIVATE(l,k,m,mm,uxdir,uydir,uu,rhocorner,rholocal,momentum_correction,cconv,invden,uvel_ref) &
!! !$OMP SHARED(f,uvel,udir,rho0,rho_relax)
!!    do k = 1,nz
!!    do l = 1,nl
!! #endif
!! 
!!       uxdir = cos(udir*pi/180.0)
!!       uydir = sin(udir*pi/180.0)
!! 
!!       uu = uvel(k)
!! 
!!       ! Floor uconv against the bulk imposed speed so a locally small
!!       ! or zero uvel(k) (e.g. near a k=1/k=nz wall-adjacent profile)
!!       ! cannot make this corner's outflow ghost effectively frozen.
!!       ! maxval(uvel) is the natural bulk reference here since the
!!       ! corner's own uconv is already the full (unprojected) speed,
!!       ! unlike the face kernels where the floor guards against an
!!       ! oblique-angle projection going to zero.
!!       uvel_ref = maxval(uvel)
!!       uconv = max(uu, uconv_min_frac*uvel_ref)
!! 
!!       cconv  = min(1.0, max(0.0, uconv))
!!       invden = 1.0/(1.0+cconv)
!! 
!!       !------------------------------------------------------------
!!       ! Fallback / mixed-corner treatment
!!       !------------------------------------------------------------
!! 
!!       f(l,0,0,k)       = 0.5*(f(l,0,1,k)     + f(l,1,0,k))
!!       f(l,0,ny+1,k)    = 0.5*(f(l,0,ny,k)    + f(l,1,ny+1,k))
!!       f(l,nx+1,0,k)    = 0.5*(f(l,nx+1,1,k)  + f(l,nx,0,k))
!!       f(l,nx+1,ny+1,k) = 0.5*(f(l,nx+1,ny,k) + f(l,nx,ny+1,k))
!! 
!!       !------------------------------------------------------------
!!       ! Quadrant 1: ux>0, uy>0  -> inflow corner (0,0), ref node (1,1,k)
!!       !------------------------------------------------------------
!!       if (uxdir > udir_tol .and. uydir > udir_tol) then
!! 
!!          if (cxs(l) >= 0 .and. cys(l) >= 0) then
!!             rholocal = 0.0
!!             do mm = 1, nl
!!                rholocal = rholocal + f(mm,1,1,k)
!!             enddo
!!             rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
!! 
!!             m = bounce(l)
!!             momentum_correction =                                     &
!!                2.0*weights(l)*rhocorner*uu                            &
!!                * (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
!!             f(l,0,0,k) = f(m,1,1,k) + momentum_correction
!!          endif
!! 
!!          f(l,nx+1,ny+1,k) =                                           &
!!             (f(l,nx+1,ny+1,k) + cconv*f(l,nx,ny,k)) * invden
!! 
!!       !------------------------------------------------------------
!!       ! Quadrant 4: ux>0, uy<0  -> inflow corner (0,ny+1), ref node (1,ny,k)
!!       !------------------------------------------------------------
!!       elseif (uxdir > udir_tol .and. uydir < -udir_tol) then
!! 
!!          if (cxs(l) >= 0 .and. cys(l) <= 0) then
!!             rholocal = 0.0
!!             do mm = 1, nl
!!                rholocal = rholocal + f(mm,1,ny,k)
!!             enddo
!!             rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
!! 
!!             m = bounce(l)
!!             momentum_correction =                                     &
!!                2.0*weights(l)*rhocorner*uu                            &
!!                * (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
!!             f(l,0,ny+1,k) = f(m,1,ny,k) + momentum_correction
!!          endif
!! 
!!          f(l,nx+1,0,k) =                                              &
!!             (f(l,nx+1,0,k) + cconv*f(l,nx,1,k)) * invden
!! 
!!       !------------------------------------------------------------
!!       ! Quadrant 2: ux<0, uy>0  -> inflow corner (nx+1,0), ref node (nx,1,k)
!!       !------------------------------------------------------------
!!       elseif (uxdir < -udir_tol .and. uydir > udir_tol) then
!! 
!!          if (cxs(l) <= 0 .and. cys(l) >= 0) then
!!             rholocal = 0.0
!!             do mm = 1, nl
!!                rholocal = rholocal + f(mm,nx,1,k)
!!             enddo
!!             rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
!! 
!!             m = bounce(l)
!!             momentum_correction =                                     &
!!                2.0*weights(l)*rhocorner*uu                            &
!!                * (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
!!             f(l,nx+1,0,k) = f(m,nx,1,k) + momentum_correction
!!          endif
!! 
!!          f(l,0,ny+1,k) =                                              &
!!             (f(l,0,ny+1,k) + cconv*f(l,1,ny,k)) * invden
!! 
!!       !------------------------------------------------------------
!!       ! Quadrant 3: ux<0, uy<0  -> inflow corner (nx+1,ny+1), ref node (nx,ny,k)
!!       !------------------------------------------------------------
!!       elseif (uxdir < -udir_tol .and. uydir < -udir_tol) then
!! 
!!          if (cxs(l) <= 0 .and. cys(l) <= 0) then
!!             rholocal = 0.0
!!             do mm = 1, nl
!!                rholocal = rholocal + f(mm,nx,ny,k)
!!             enddo
!!             rhocorner = rho_relax*rholocal + (1.0-rho_relax)*rho0
!! 
!!             m = bounce(l)
!!             momentum_correction =                                     &
!!                2.0*weights(l)*rhocorner*uu                            &
!!                * (real(cxs(l))*uxdir + real(cys(l))*uydir)/cs2
!!             f(l,nx+1,ny+1,k) = f(m,nx,ny,k) + momentum_correction
!!          endif
!! 
!!          f(l,0,0,k) =                                                 &
!!             (f(l,0,0,k) + cconv*f(l,1,1,k)) * invden
!! 
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
