module m_boundary_j_k_edges
contains
subroutine boundary_j_k_edges(f1,f2,opt)
   use mod_dimensions
   use mod_D3Q27setup, only : nl, cxs, cys, czs
   implicit none
   real, intent(inout) :: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, value :: opt   ! -1 = no-slip, +1 = free-slip
   integer :: i,j,k,l,m

#ifdef _CUDA
   attributes(device) :: f1, f2
#endif

!======================================================================
! j=1, k=1 edge → ghost (j=0, k=0)
! normals: y and z → reverse cy, cz
! tangential: x stays same (always)
!====================================================================

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) < 0 .and. czs(l) < 0) then
            do m=1,nl
               if ( cxs(m) == opt*cxs(l) .and. &
                    cys(m) ==    -cys(l) .and. &
                    czs(m) ==    -czs(l) ) then
                   f1(m,i,0,0) = f2(l,i,0,0)
                   exit
               endif
            enddo

         endif
      enddo
   enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) < 0 .and. czs(l) < 0) then
            f1(l,i,0,0) = f1(l, i, 1, 1)   ! interior: only y,z shift
         endif
      enddo
   enddo


!====================================================================
! j=ny, k=1 → ghost (ny+1, 0)
!====================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) > 0 .and. czs(l) < 0) then
            do m=1,nl
               if ( cxs(m)== opt*cxs(l) .and. &
                    cys(m)==    -cys(l) .and. &
                    czs(m)==    -czs(l)) then
                  f1(m,i,ny+1,0) = f2(l,i,ny+1,0)
                  exit
               endif
            enddo
         endif
      enddo
   enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) > 0 .and. czs(l) < 0) then
            f1(l,i,ny+1,0) = f1(l, i, ny, 1)
         endif
      enddo
   enddo


!====================================================================
! j=1, k=nz → ghost (0, nz+1)
!====================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) < 0 .and. czs(l) > 0) then
            do m=1,nl
               if (cxs(m)== opt*cxs(l) .and. &
                   cys(m)==    -cys(l) .and. &
                   czs(m)==    -czs(l)) then
                  f1(m,i,0,nz+1)=f2(l,i,0,nz+1)
                  exit
               endif
            enddo
         endif
      enddo
   enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) < 0 .and. czs(l) > 0) then
            f1(l,i,0,nz+1) = f1(l, i, 1, nz)
         endif
      enddo
   enddo


!====================================================================
! j=ny, k=nz → ghost (ny+1, nz+1)
!====================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) > 0 .and. czs(l) > 0) then
            do m=1,nl
               if (cxs(m)== cxs(l) .and. &
                   cys(m)==-cys(l) .and. &
                   czs(m)==-czs(l)) then
                  f1(m,i,ny+1,nz+1)=f2(l,i,ny+1,nz+1)
                  exit
               endif
            enddo
         endif
      enddo
   enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) > 0 .and. czs(l) > 0) then
            f1(l,i,ny+1,nz+1) = f1(l, i, ny, nz)
         endif
      enddo
   enddo

end subroutine
end module

