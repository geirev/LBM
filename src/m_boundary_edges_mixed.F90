module m_boundary_edges_mixed
contains

subroutine boundary_edges_mixed(f,ibnd,jbnd,kbnd,j0_is_phys,jN_is_phys)

!-----------------------------------------------------------------------
! Completion of mixed open/closed boundary edges.
!
! This routine is consistent with:
!
!   boundary_inflow_i_kernel
!   boundary_inflow_j_kernel
!   boundary_inflow_edges_ij_kernel
!
! The general rule for an edge where one direction is OPEN and the
! other direction is CLOSED is:
!
!   1. The closed-boundary routine first constructs the corresponding
!      wall ghost plane.
!
!   2. The mixed edge is then obtained by zero-normal-gradient
!      extrapolation in the OPEN direction.
!
! Thus no history-dependent convective term is used anywhere.
!
! Examples:
!
!   i open, j closed:
!
!      f(:,0,0,k) = f(:,1,0,k)
!
!   The j-wall value at j=0 is retained, while the value is simply
!   extrapolated in the open i-direction.
!
!
!   j open, i closed:
!
!      f(:,0,0,k) = f(:,0,1,k)
!
!   The i-wall value at i=0 is retained, while the value is simply
!   extrapolated in the open j-direction.
!
!
! Likewise for mixed i-k and j-k edges.
!
!
! IMPORTANT:
!
! This routine assumes that the face boundary conditions have already
! been applied. In particular, the CLOSED boundary ghost planes must
! already contain their final wall values before this routine is called.
!
!
! MPI:
!
!   j0_is_phys = .true. only when local j=0 is a physical boundary.
!   jN_is_phys = .true. only when local j=ny+1 is a physical boundary.
!
! Any edge touching a non-physical j halo is left untouched. The MPI
! halo exchange is responsible for those values.
!
!
! This routine does NOT handle the case i-open + j-open. Those four
! corners are handled by boundary_inflow_edges_ij_kernel.
!
!-----------------------------------------------------------------------

   use mod_dimensions
   use mod_D3Q27setup, only : nl

   implicit none

   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)

   integer, intent(in) :: ibnd
   integer, intent(in) :: jbnd
   integer, intent(in) :: kbnd

   logical, intent(in) :: j0_is_phys
   logical, intent(in) :: jN_is_phys

#ifdef _CUDA
   attributes(device) :: f
#endif

   integer :: l,i,j,k


!=======================================================================
! i OPEN, j CLOSED
!
! The j-wall ghost values have already been constructed.
! Extrapolate those values in the open i-direction.
!
! Lower j boundary:
!
!     (1,0)  --> (0,0)
!     (nx,0) --> (nx+1,0)
!
! Upper j boundary:
!
!     (1,ny+1)  --> (0,ny+1)
!     (nx,ny+1) --> (nx+1,ny+1)
!
! No udir-dependent weighting is needed: the j-boundary is closed and
! the i-boundary uses zero-gradient whenever it is not imposing inflow.
!=======================================================================

   if (ibnd == 1 .and. jbnd > 10) then

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k = 1,nz
      do l = 1,nl

         if (j0_is_phys) then

            f(l,0,0,k)    = f(l,1,0,k)
            f(l,nx+1,0,k) = f(l,nx,0,k)

         endif


         if (jN_is_phys) then

            f(l,0,ny+1,k)    = f(l,1,ny+1,k)
            f(l,nx+1,ny+1,k) = f(l,nx,ny+1,k)

         endif

      enddo
      enddo

   endif


!=======================================================================
! j OPEN, i CLOSED
!
! The i-wall ghost values have already been constructed.
! Extrapolate those values in the open j-direction.
!
! Lower j boundary:
!
!     (0,1)    --> (0,0)
!     (nx+1,1) --> (nx+1,0)
!
! Upper j boundary:
!
!     (0,ny)    --> (0,ny+1)
!     (nx+1,ny) --> (nx+1,ny+1)
!
! Only genuine physical j boundaries are modified.
!=======================================================================

   if (jbnd == 1 .and. ibnd > 10) then

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k = 1,nz
      do l = 1,nl

         if (j0_is_phys) then

            f(l,0,0,k)    = f(l,0,1,k)
            f(l,nx+1,0,k) = f(l,nx+1,1,k)

         endif


         if (jN_is_phys) then

            f(l,0,ny+1,k)    = f(l,0,ny,k)
            f(l,nx+1,ny+1,k) = f(l,nx+1,ny,k)

         endif

      enddo
      enddo

   endif


!=======================================================================
! i OPEN, k CLOSED
!
! The k-wall ghost values have already been constructed.
! Extrapolate those values in the open i-direction.
!
! Bottom k boundary:
!
!     (1,j,0)  --> (0,j,0)
!     (nx,j,0) --> (nx+1,j,0)
!
! Top k boundary:
!
!     (1,j,nz+1)  --> (0,j,nz+1)
!     (nx,j,nz+1) --> (nx+1,j,nz+1)
!
! Only j=1..ny is handled here. The true i-j-k corners are deliberately
! excluded.
!=======================================================================

   if (ibnd == 1 .and. kbnd > 10) then

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do j = 1,ny
      do l = 1,nl

         f(l,0,j,0)    = f(l,1,j,0)
         f(l,nx+1,j,0) = f(l,nx,j,0)

         f(l,0,j,nz+1)    = f(l,1,j,nz+1)
         f(l,nx+1,j,nz+1) = f(l,nx,j,nz+1)

      enddo
      enddo

   endif


!=======================================================================
! j OPEN, k CLOSED
!
! The k-wall ghost values have already been constructed.
! Extrapolate those values in the open j-direction.
!
! Lower j boundary:
!
!     (i,1,0)    --> (i,0,0)
!     (i,1,nz+1) --> (i,0,nz+1)
!
! Upper j boundary:
!
!     (i,ny,0)    --> (i,ny+1,0)
!     (i,ny,nz+1) --> (i,ny+1,nz+1)
!
! Only genuine physical j boundaries are modified.
!
! i=1..nx only; the true i-j-k corners are deliberately excluded.
!=======================================================================

   if (jbnd == 1 .and. kbnd > 10) then

#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do i = 1,nx
      do l = 1,nl

         if (j0_is_phys) then

            f(l,i,0,0)    = f(l,i,1,0)
            f(l,i,0,nz+1) = f(l,i,1,nz+1)

         endif


         if (jN_is_phys) then

            f(l,i,ny+1,0)    = f(l,i,ny,0)
            f(l,i,ny+1,nz+1) = f(l,i,ny,nz+1)

         endif

      enddo
      enddo

   endif


end subroutine boundary_edges_mixed

end module m_boundary_edges_mixed
!!!module m_boundary_edges_mixed
!!!contains
!!!subroutine boundary_edges_mixed(f,ibnd,jbnd,kbnd,j0_is_phys,jN_is_phys)
!!!! Mixed open-closed edge completion.
!!!!
!!!!  - j0_is_phys / jN_is_phys: true only if the local j=0 / j=ny+1 plane
!!!!    is a genuine physical domain boundary on this MPI rank (false at
!!!!    an inter-rank halo interface, where the halo exchange - not this
!!!!    routine - is responsible for that plane's content). Only the two
!!!!    branches below that fire on jbnd==1 need this; the ibnd==1 branches
!!!!    are unaffected since i is never MPI-decomposed here, and the
!!!!    jbnd>10 (closed) branches assume j is not decomposed under a
!!!!    closed boundary - if you ever run MPI-decomposed j together with
!!!!    a closed jbnd, that combination needs separate treatment not
!!!!    covered here.
!!!   use mod_dimensions
!!!   use mod_D3Q27setup, only : nl
!!!   implicit none
!!!   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
!!!   integer, intent(in) :: ibnd,jbnd,kbnd
!!!   logical, intent(in) :: j0_is_phys, jN_is_phys
!!!#ifdef _CUDA
!!!   attributes(device) :: f
!!!#endif
!!!   integer :: l,i,j,k
!!!
!!!   ! i open, j closed: extrapolate in i on the completed j-wall ghosts.
!!!   ! Not j0/jN-gated: jbnd>10 here, so j is a closed wall, not an
!!!   ! open-boundary plane subject to the halo-interface ambiguity.
!!!   if (ibnd==1 .and. jbnd>10) then
!!!#ifdef _CUDA
!!!!$cuf kernel do(2) <<<*,*>>>
!!!#endif
!!!      do k=1,nz
!!!      do l=1,nl
!!!         f(l,0,0,k)       = f(l,1,0,k)
!!!         f(l,nx+1,0,k)    = f(l,nx,0,k)
!!!         f(l,0,ny+1,k)    = f(l,1,ny+1,k)
!!!         f(l,nx+1,ny+1,k) = f(l,nx,ny+1,k)
!!!      enddo
!!!      enddo
!!!   endif
!!!
!!!   ! j open, i closed: extrapolate in j on the completed i-wall ghosts.
!!!   if (jbnd==1 .and. ibnd>10) then
!!!#ifdef _CUDA
!!!!$cuf kernel do(2) <<<*,*>>>
!!!#endif
!!!      do k=1,nz
!!!      do l=1,nl
!!!         if (j0_is_phys) then
!!!            f(l,0,0,k)    = f(l,0,1,k)
!!!            f(l,nx+1,0,k) = f(l,nx+1,1,k)
!!!         endif
!!!         if (jN_is_phys) then
!!!            f(l,0,ny+1,k)    = f(l,0,ny,k)
!!!            f(l,nx+1,ny+1,k) = f(l,nx+1,ny,k)
!!!         endif
!!!      enddo
!!!      enddo
!!!   endif
!!!
!!!   ! i open, k closed.
!!!   ! Not j0/jN-gated: this branch never touches a j=0/ny+1 plane.
!!!   if (ibnd==1 .and. kbnd>10) then
!!!#ifdef _CUDA
!!!!$cuf kernel do(2) <<<*,*>>>
!!!#endif
!!!      do j=1,ny
!!!      do l=1,nl
!!!         f(l,0,j,0)       = f(l,1,j,0)
!!!         f(l,nx+1,j,0)    = f(l,nx,j,0)
!!!         f(l,0,j,nz+1)    = f(l,1,j,nz+1)
!!!         f(l,nx+1,j,nz+1) = f(l,nx,j,nz+1)
!!!      enddo
!!!      enddo
!!!   endif
!!!
!!!   ! j open, k closed.
!!!   if (jbnd==1 .and. kbnd>10) then
!!!#ifdef _CUDA
!!!!$cuf kernel do(2) <<<*,*>>>
!!!#endif
!!!      do i=1,nx
!!!      do l=1,nl
!!!         if (j0_is_phys) then
!!!            f(l,i,0,0)    = f(l,i,1,0)
!!!            f(l,i,0,nz+1) = f(l,i,1,nz+1)
!!!         endif
!!!         if (jN_is_phys) then
!!!            f(l,i,ny+1,0)    = f(l,i,ny,0)
!!!            f(l,i,ny+1,nz+1) = f(l,i,ny,nz+1)
!!!         endif
!!!      enddo
!!!      enddo
!!!   endif
!!!end subroutine
!!!end module
