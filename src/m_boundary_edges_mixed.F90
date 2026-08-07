module m_boundary_edges_mixed
contains
subroutine boundary_edges_mixed(f,ibnd,jbnd,kbnd,j0_is_phys,jN_is_phys)
! Mixed open-closed edge completion.
!
!  - j0_is_phys / jN_is_phys: true only if the local j=0 / j=ny+1 plane
!    is a genuine physical domain boundary on this MPI rank (false at
!    an inter-rank halo interface, where the halo exchange - not this
!    routine - is responsible for that plane's content). Only the two
!    branches below that fire on jbnd==1 need this; the ibnd==1 branches
!    are unaffected since i is never MPI-decomposed here, and the
!    jbnd>10 (closed) branches assume j is not decomposed under a
!    closed boundary - if you ever run MPI-decomposed j together with
!    a closed jbnd, that combination needs separate treatment not
!    covered here.
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, intent(in) :: ibnd,jbnd,kbnd
   logical, intent(in) :: j0_is_phys, jN_is_phys
#ifdef _CUDA
   attributes(device) :: f
#endif
   integer :: l,i,j,k

   ! i open, j closed: extrapolate in i on the completed j-wall ghosts.
   ! Not j0/jN-gated: jbnd>10 here, so j is a closed wall, not an
   ! open-boundary plane subject to the halo-interface ambiguity.
   if (ibnd==1 .and. jbnd>10) then
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=1,nz
      do l=1,nl
         f(l,0,0,k)       = f(l,1,0,k)
         f(l,nx+1,0,k)    = f(l,nx,0,k)
         f(l,0,ny+1,k)    = f(l,1,ny+1,k)
         f(l,nx+1,ny+1,k) = f(l,nx,ny+1,k)
      enddo
      enddo
   endif

   ! j open, i closed: extrapolate in j on the completed i-wall ghosts.
   if (jbnd==1 .and. ibnd>10) then
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=1,nz
      do l=1,nl
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

   ! i open, k closed.
   ! Not j0/jN-gated: this branch never touches a j=0/ny+1 plane.
   if (ibnd==1 .and. kbnd>10) then
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do j=1,ny
      do l=1,nl
         f(l,0,j,0)       = f(l,1,j,0)
         f(l,nx+1,j,0)    = f(l,nx,j,0)
         f(l,0,j,nz+1)    = f(l,1,j,nz+1)
         f(l,nx+1,j,nz+1) = f(l,nx,j,nz+1)
      enddo
      enddo
   endif

   ! j open, k closed.
   if (jbnd==1 .and. kbnd>10) then
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do i=1,nx
      do l=1,nl
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
end subroutine
end module
