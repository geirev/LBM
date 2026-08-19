module m_boundary_inflow_j_kernel
contains

#ifdef _CUDA
attributes(global) &
#endif
subroutine boundary_inflow_j_kernel(f,uvel,udir,rho0,rho_relax, &
                                    inv1cs2,inv2cs4,inv6cs6,ibgk, &
                                    j0_is_phys,jN_is_phys)

!-----------------------------------------------------------------------
! General inflow / outflow boundary condition in j direction.
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
! non-equilibrium extrapolation:
!
!       fghost =
!
!           feq(rho0,u_local)
!
!         + [f_int - feq(rho_local,u_local)]
!
!
! WIDE COMPLEMENTARY TRANSITION
! =============================
!
! The boundary transitions directly and smoothly between
!
!       inflow <-> NEE outflow
!
! without an intermediate zero-gradient state for the populations
! entering from the ghost plane.
!
! At each physical j face:
!
!       win + wout = 1
!
! for all flow directions.
!
! At uydir = 0:
!
!       win  = 0.5
!       wout = 0.5
!
! The transition uses a quintic smootherstep:
!
!       S(s) = 6*s^5 - 15*s^4 + 10*s^3
!
! with uydir mapped across +/- blend_width.
!
! Only populations capable of streaming FROM the ghost plane INTO the
! computational domain are reconstructed:
!
!       j=0    : cy > 0
!       j=ny+1 : cy < 0
!
! All other populations remain zero-gradient copies of the adjacent
! interior populations.
!
! j0_is_phys / jN_is_phys prevent modification of internal MPI
! j interfaces.
!-----------------------------------------------------------------------

#ifdef _CUDA
   use cudafor
#endif

   use mod_dimensions
   use mod_D3Q27setup, only : nl,cxs,cys,czs,bounce

   use m_fequil_boundary, only : &
        fequil_value,fequil_difference

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

   ! Use the same width as the successful i-direction test.
   real, parameter :: blend_width = 0.30


!-----------------------------------------------------------------------
! Prescribed velocity
!-----------------------------------------------------------------------

   real :: uxdir,uydir
   real :: uu
   real :: ux,uy


!-----------------------------------------------------------------------
! Local density
!-----------------------------------------------------------------------

   real :: rholocal0
   real :: rholocalN

   real :: rholoc0
   real :: rholocN


!-----------------------------------------------------------------------
! Local momentum / velocity for NEE outflow
!-----------------------------------------------------------------------

   real :: momx0,momy0,momz0
   real :: momxN,momyN,momzN

   real :: uloc0,vloc0,wloc0
   real :: ulocN,vlocN,wlocN


!-----------------------------------------------------------------------
! Complementary transition weights
!-----------------------------------------------------------------------

   real :: x
   real :: s

   real :: win0
   real :: winN

   real :: wout0
   real :: woutN


!-----------------------------------------------------------------------
! Distribution values
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
!$OMP&         x,s,win0,winN,wout0,woutN,                          &
!$OMP&         fint0,fintN,fin0,finN,fout0,foutN,                  &
!$OMP&         feq_int,feq_out,dfeq)                               &
!$OMP& SHARED(f,uvel,udir,rho0,rho_relax,                          &
!$OMP&        inv1cs2,inv2cs4,inv6cs6,ibgk,                       &
!$OMP&        j0_is_phys,jN_is_phys)

   do k = 1,nz
   do i = 1,nx

#endif


!=======================================================================
! Prescribed undisturbed velocity
!=======================================================================

      uxdir = cos(udir*pi/180.0)
      uydir = sin(udir*pi/180.0)

      uu = uvel(k)

      ux = uu*uxdir
      uy = uu*uydir


!=======================================================================
! Complementary smooth inflow/outflow weights.
!
! Map uydir to x in [-1,1]:
!
!       uydir <= -blend_width : x = -1
!       uydir = 0             : x =  0
!       uydir >= +blend_width : x = +1
!
! Then map to s in [0,1] and apply quintic smootherstep.
!
!
! LOWER FACE j=0:
!
!       uy << 0 : outflow
!       uy = 0  : 50% inflow / 50% outflow
!       uy >> 0 : inflow
!
!
! UPPER FACE j=ny+1:
!
!       uy << 0 : inflow
!       uy = 0  : 50% inflow / 50% outflow
!       uy >> 0 : outflow
!=======================================================================

      x = min(1.0,max(-1.0,uydir/blend_width))

      s = 0.5*(x+1.0)

      s = s*s*s*(10.0 + s*(-15.0 + 6.0*s))


!-----------------------------------------------------------------------
! Lower face j=0
!-----------------------------------------------------------------------

      win0  = s
      wout0 = 1.0-s


!-----------------------------------------------------------------------
! Upper face j=ny+1
!-----------------------------------------------------------------------

      winN  = 1.0-s
      woutN = s


!=======================================================================
! Local density and momentum at adjacent interior nodes
!=======================================================================

      rholocal0 = 0.0

      momx0 = 0.0
      momy0 = 0.0
      momz0 = 0.0


      rholocalN = 0.0

      momxN = 0.0
      momyN = 0.0
      momzN = 0.0


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
! Local velocity used for NEE outflow
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
! Default for populations not entering from the ghost plane:
! zero-gradient.
!-----------------------------------------------------------------------

            f(l,i,0,k) = fint0


!-----------------------------------------------------------------------
! Only cy>0 populations can enter through j=0.
!-----------------------------------------------------------------------

            if (cys(l) > 0) then


!=======================================================================
! LOWER INFLOW CANDIDATE
!=======================================================================

               m = bounce(l)

               dfeq = fequil_difference( &
                          l,rholoc0,ux,uy,0.0, &
                          inv1cs2,inv6cs6,ibgk)

               fin0 = f(m,i,1,k) + dfeq


!=======================================================================
! LOWER OUTFLOW CANDIDATE
!=======================================================================

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


!=======================================================================
! Complementary blend
!=======================================================================

               f(l,i,0,k) = &
                    win0 *fin0 + &
                    wout0*fout0

            endif

         endif


!#######################################################################
! UPPER FACE: j=ny+1
!#######################################################################

         if (jN_is_phys) then

            fintN = f(l,i,ny,k)


!-----------------------------------------------------------------------
! Default for populations not entering from the ghost plane:
! zero-gradient.
!-----------------------------------------------------------------------

            f(l,i,ny+1,k) = fintN


!-----------------------------------------------------------------------
! Only cy<0 populations can enter through j=ny+1.
!-----------------------------------------------------------------------

            if (cys(l) < 0) then


!=======================================================================
! UPPER INFLOW CANDIDATE
!=======================================================================

               m = bounce(l)

               dfeq = fequil_difference( &
                          l,rholocN,ux,uy,0.0, &
                          inv1cs2,inv6cs6,ibgk)

               finN = f(m,i,ny,k) + dfeq


!=======================================================================
! UPPER OUTFLOW CANDIDATE
!=======================================================================

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


!=======================================================================
! Complementary blend
!=======================================================================

               f(l,i,ny+1,k) = &
                    winN *finN + &
                    woutN*foutN

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
