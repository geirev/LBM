module m_boundary_i_j_edges
contains
subroutine boundary_i_j_edges(f1,f2,opt)
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
! i=1, j=1 → ghost (0,0)
! normals: x, y  → reverse cx, cy
! tangential: z  → cz reversed only for no-slip
!====================================================================

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do k=1,nz
      do l=1,nl
         if (cxs(l) < 0 .and. cys(l) < 0) then
            do m=1,nl
               if ( cxs(m) == -cxs(l) .and. &
                    cys(m) == -cys(l) .and. &
                    czs(m) == opt*czs(l) ) then
                  f1(m,0,0,k) = f2(l,0,0,k)
                  exit
               endif
            enddo
         endif
      enddo
   enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do k=1,nz
      do l=1,nl
         if (cxs(l) < 0 .and. cys(l) < 0) then
            f1(l,0,0,k) = f1(l, 1, 1, k)
         endif
      enddo
   enddo


!====================================================================
! i=nx, j=1 → ghost (nx+1,0)
!====================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do k=1,nz
      do l=1,nl
         if (cxs(l) > 0 .and. cys(l) < 0) then

            do m=1,nl
               if ( cxs(m)==-cxs(l) .and. &
                    cys(m)==-cys(l) .and. &
                    czs(m)== opt*czs(l)) then
                  f1(m,nx+1,0,k) = f2(l,nx+1,0,k)
                  exit
               endif
            enddo
         endif
      enddo
   enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do k=1,nz
      do l=1,nl
         if (cxs(l) > 0 .and. cys(l) < 0) then
            f1(l,nx+1,0,k) = f1(l, nx, 1, k)
         endif
      enddo
   enddo


!====================================================================
! i=1, j=ny → ghost (0,ny+1)
!====================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do k=1,nz
      do l=1,nl
         if (cxs(l)<0 .and. cys(l)>0) then
            do m=1,nl
               if (cxs(m)==-cxs(l) .and. cys(m)==-cys(l) .and. czs(m)==opt*czs(l)) then
                  f1(m,0,ny+1,k) = f2(l,0,ny+1,k)
                  exit
               endif
            enddo
         endif
      enddo
   enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do k=1,nz
      do l=1,nl
         if (cxs(l)<0 .and. cys(l)>0) then
            f1(l,0,ny+1,k) = f1(l, 1, ny, k)
         endif
      enddo
   enddo


!======================================================================
! i=nx, j=ny → ghost (nx+1,ny+1)
!======================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do k=1,nz
      do l=1,nl
         if (cxs(l)>0 .and. cys(l)>0) then
            do m=1,nl
               if (cxs(m)==-cxs(l) .and. cys(m)==-cys(l) .and. czs(m)==opt*czs(l)) then
                  f1(m,nx+1,ny+1,k) = f2(l,nx+1,ny+1,k)
                  exit
               endif
            enddo
         endif
      enddo
   enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do k=1,nz
      do l=1,nl
         if (cxs(l)>0 .and. cys(l)>0) then
            f1(l,nx+1,ny+1,k) = f1(l, nx, ny, k)
         endif
      enddo
   enddo

end subroutine
end module

