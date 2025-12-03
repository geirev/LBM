module m_boundary_i_k_edges
contains
subroutine boundary_i_k_edges(f1,f2,opt)
   use mod_dimensions
   use mod_D3Q27setup, only : nl, cxs, cys, czs
   implicit none
   real, intent(inout) :: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, value :: opt    ! -1 no-slip, +1 free-slip
   integer :: i,j,k,l,m

#ifdef _CUDA
   attributes(device) :: f1, f2
#endif

!======================================================================
! i=1, k=1  → ghost (0, j, 0)
! normals: x, z → reverse cx, cz
! y is tangential → cy preserved for free-slip, reversed for no-slip
!======================================================================

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
do j=1,ny
   do l=1,nl
      if (cxs(l)<0 .and. czs(l)<0) then

         do m=1,nl
            if ( cxs(m)==-cxs(l) .and. &
                 cys(m)== opt*cys(l) .and. &
                 czs(m)==-czs(l) ) then
               f1(m,0,j,0)=f2(l,0,j,0)
               exit
            endif
         enddo

      endif
   enddo
enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
do j=1,ny
   do l=1,nl
      if (cxs(l)<0 .and. czs(l)<0) then
         f1(l,0,j,0) = f1(l, 1, j, 1)
      endif
   enddo
enddo


!======================================================================
! i=nx, k=1 → (nx+1,j,0)
!======================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
do j=1,ny
   do l=1,nl
      if (cxs(l)>0 .and. czs(l)<0) then

         do m=1,nl
            if (cxs(m)==-cxs(l) .and. cys(m)==opt*cys(l) .and. czs(m)==-czs(l)) then
               f1(m,nx+1,j,0)=f2(l,nx+1,j,0)
               exit
            endif
         enddo

      endif
   enddo
enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
do j=1,ny
   do l=1,nl
      if (cxs(l)>0 .and. czs(l)<0) then
         f1(l,nx+1,j,0)=f1(l, nx, j, 1)
      endif
   enddo
enddo


!======================================================================
! i=1, k=nz  → (0,j,nz+1)
!======================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
do j=1,ny
   do l=1,nl
      if (cxs(l)<0 .and. czs(l)>0) then

         do m=1,nl
            if (cxs(m)==-cxs(l) .and. cys(m)==opt*cys(l) .and. czs(m)==-czs(l)) then
               f1(m,0,j,nz+1)=f2(l,0,j,nz+1)
               exit
            endif
         enddo

      endif
   enddo
enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
do j=1,ny
   do l=1,nl
      if (cxs(l)<0 .and. czs(l)>0) then
         f1(l,0,j,nz+1)=f1(l,1,j,nz)
      endif
   enddo
enddo


!======================================================================
! i=nx, k=nz → (nx+1,j,nz+1)
!======================================================================
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
do j=1,ny
   do l=1,nl
      if (cxs(l)>0 .and. czs(l)>0) then

         do m=1,nl
            if (cxs(m)==-cxs(l) .and. cys(m)==opt*cys(l) .and. czs(m)==-czs(l)) then
               f1(m,nx+1,j,nz+1)=f2(l,nx+1,j,nz+1)
               exit
            endif
         enddo

      endif
   enddo
enddo

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
do j=1,ny
   do l=1,nl
      if (cxs(l)>0 .and. czs(l)>0) then
         f1(l,nx+1,j,nz+1)=f1(l,nx,j,nz)
      endif
   enddo
enddo

end subroutine
end module

