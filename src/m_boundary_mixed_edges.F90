module m_boundary_mixed_edges
contains
subroutine boundary_mixed_edges(f,ibnd,jbnd,kbnd)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, intent(in) :: ibnd,jbnd,kbnd
#ifdef _CUDA
   attributes(device) :: f
#endif
   integer :: l,i,j,k

   ! i open, j closed: extrapolate in i on the completed j-wall ghosts.
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
         f(l,0,0,k)       = f(l,0,1,k)
         f(l,0,ny+1,k)    = f(l,0,ny,k)
         f(l,nx+1,0,k)    = f(l,nx+1,1,k)
         f(l,nx+1,ny+1,k) = f(l,nx+1,ny,k)
      enddo
      enddo
   endif

   ! i open, k closed.
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
         f(l,i,0,0)       = f(l,i,1,0)
         f(l,i,ny+1,0)    = f(l,i,ny,0)
         f(l,i,0,nz+1)    = f(l,i,1,nz+1)
         f(l,i,ny+1,nz+1) = f(l,i,ny,nz+1)
      enddo
      enddo
   endif
end subroutine boundary_mixed_edges
end module m_boundary_mixed_edges
