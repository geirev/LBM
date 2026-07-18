module m_boundary_nongeneric_corners
contains
subroutine boundary_nongeneric_corners(f,ibnd,jbnd,kbnd)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, intent(in) :: ibnd,jbnd,kbnd
#ifdef _CUDA
   attributes(device) :: f
#endif
   integer :: l

   ! Periodic intersections are finalized later by periodic sweeps.
   if (ibnd==0 .or. jbnd==0 .or. kbnd==0) return

   ! The only currently supported nonperiodic three-dimensional case
   ! involving open boundaries has kbnd closed and ibnd/jbnd open or closed.
   if (kbnd>10) then
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do l=1,nl
         f(l,0,0,0) = 0.5*(f(l,1,0,0)+f(l,0,1,0))
         f(l,nx+1,0,0) = 0.5*(f(l,nx,0,0)+f(l,nx+1,1,0))
         f(l,0,ny+1,0) = 0.5*(f(l,1,ny+1,0)+f(l,0,ny,0))
         f(l,nx+1,ny+1,0) = 0.5*(f(l,nx,ny+1,0)+f(l,nx+1,ny,0))

         f(l,0,0,nz+1) = 0.5*(f(l,1,0,nz+1)+f(l,0,1,nz+1))
         f(l,nx+1,0,nz+1) = 0.5*(f(l,nx,0,nz+1)+f(l,nx+1,1,nz+1))
         f(l,0,ny+1,nz+1) = 0.5*(f(l,1,ny+1,nz+1)+f(l,0,ny,nz+1))
         f(l,nx+1,ny+1,nz+1) = 0.5*(f(l,nx,ny+1,nz+1)+f(l,nx+1,ny,nz+1))
      enddo
   endif
end subroutine boundary_nongeneric_corners
end module m_boundary_nongeneric_corners

