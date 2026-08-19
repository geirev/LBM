module m_boundary_inflow_edges_ij_kernel
contains

#ifdef _CUDA
attributes(global) &
#endif
subroutine boundary_inflow_edges_ij_kernel( &
   f,uvel,udir,rho0,rho_relax,inv1cs2,inv6cs6,ibgk, &
   j0_is_phys,jN_is_phys)

!-----------------------------------------------------------------------
! General treatment of the four true i-j edges
! (intersection of two open side boundaries), for k = 1..nz.
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
! ===================
!
! Each edge is first constructed from the two adjacent completed face
! ghost distributions:
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
! Since the face ghost values already contain the complementary
! inflow/NEE-outflow blending, this baseline automatically inherits
! the corresponding face treatment.
!
!
! TRUE INFLOW-INFLOW POPULATIONS
! ==============================
!
! Only populations crossing BOTH boundaries are specially
! reconstructed:
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
! where feq(l)-feq(bounce(l)) uses exactly the same equilibrium order
! as postcoll_kernel.
!
!
! COMPLEMENTARY FACE WEIGHTS
! ==========================
!
! The i and j face-inflow weights are constructed using exactly the
! same wide complementary quintic smootherstep as the face kernels.
!
! For i:
!
!    win_i0 + win_iN = 1
!
! and for j:
!
!    win_j0 + win_jN = 1
!
! throughout the transition.
!
! At ux = 0:
!
!    win_i0 = win_iN = 0.5
!
! At uy = 0:
!
!    win_j0 = win_jN = 0.5
!
! The true two-face inflow reconstruction is then weighted by
!
!    win00 = win_i0 * win_j0
!    win0N = win_i0 * win_jN
!    winN0 = win_iN * win_j0
!    winNN = win_iN * win_jN
!
! This makes the edge treatment consistent with the wide smooth
! inflow/outflow role transition of the face kernels.
!
!
! MPI
! ===
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

   use m_fequil_boundary, only : fequil_difference

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
   real, value :: inv6cs6

   integer, value :: ibgk

   logical, value :: j0_is_phys
   logical, value :: jN_is_phys


#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uvel
#endif


!=======================================================================
! Local variables
!=======================================================================

   integer :: k,l,m,mm

   real, parameter :: pi = acos(-1.0)

   ! Must match boundary_inflow_i_kernel and
   ! boundary_inflow_j_kernel.
   real, parameter :: blend_width = 0.30

   real :: uxdir,uydir
   real :: ax,ay,aden

   real :: uu
   real :: ux,uy

   real :: sx,sy

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
!$OMP& PRIVATE(l,k,m,mm,uxdir,uydir,ax,ay,aden,uu,ux,uy,sx,sy,      &
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

      ! uxdir^2 + uydir^2 = 1, so this is always positive.
      aden = ax + ay


!=======================================================================
! Complementary i-face weights
!
! Exactly the same mapping as boundary_inflow_i_kernel.
!
! First map:
!
!      uxdir/blend_width -> [-1,1]
!
! then to [0,1], followed by quintic smootherstep:
!
!      S(s) = 6*s^5 - 15*s^4 + 10*s^3
!
! At uxdir = 0:
!
!      win_i0 = 0.5
!      win_iN = 0.5
!=======================================================================

      sx = min(1.0,max(-1.0,uxdir/blend_width))

      sx = 0.5*(sx+1.0)

      sx = sx*sx*sx*(10.0 + sx*(-15.0 + 6.0*sx))

      win_i0 = sx
      win_iN = 1.0-sx


!=======================================================================
! Complementary j-face weights
!
! Exactly the same mapping as boundary_inflow_j_kernel.
!
! At uydir = 0:
!
!      win_j0 = 0.5
!      win_jN = 0.5
!=======================================================================

      sy = min(1.0,max(-1.0,uydir/blend_width))

      sy = 0.5*(sy+1.0)

      sy = sy*sy*sy*(10.0 + sy*(-15.0 + 6.0*sy))

      win_j0 = sy
      win_jN = 1.0-sy


!=======================================================================
! True two-face inflow weights
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
