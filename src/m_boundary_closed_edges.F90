module m_boundary_closed_edges
contains
subroutine boundary_closed_edges(f1,f2,opt)
   use mod_dimensions
   use mod_D3Q27setup, only : cxs,cys,czs,nl
   implicit none
   real, intent(inout):: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout):: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, value     :: opt    ! -1 for noslip, 1 for freeslip
#ifdef _CUDA
   attributes(device) :: f1
   attributes(device) :: f2
#endif
   integer i,j,k,l,m

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! === Edge along j=1 and k=1 ===
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if ((cys(l) < 0) .and. (czs(l) < 0)) then
            do m=1,nl
               if (cxs(m) == (opt*cxs(l)) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
                  f1(m,i,0,0)=f2(l,i,0,0)
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
         if ((cys(l) < 0) .and. (czs(l) < 0)) then
            j=0; k=0
            f1(l, i, j, k)=f1(l, i-cxs(l), j-cys(l), k-czs(l))
         endif
      enddo
   enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! === Edge along j=ny and k=1 ===
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) > 0 .and. czs(l) < 0) then
            do m=1,nl
               if (cxs(m) == (opt*cxs(l)) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
                  f1(m,i,ny+1,0)=f2(l,i,ny+1,0)
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
            j=ny+1; k=0
            f1(l, i, j, k)=f1(l, i-cxs(l), j-cys(l), k-czs(l))
         endif
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! === Edge along j=1 and k=nz ===
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) < 0 .and. czs(l) > 0) then
            do m=1,nl
               if (cxs(m) == (opt*cxs(l)) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
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
            j=0; k=nz+1
            f1(l, i, j, k)=f1(l, i-cxs(l), j-cys(l), k-czs(l))
         endif
      enddo
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! === Edge along j=ny and k=nz ===
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
   do i=1,nx
      do l=1,nl
         if (cys(l) > 0 .and. czs(l) > 0) then
            do m=1,nl
               if (cxs(m) == (opt*cxs(l)) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
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
            j=ny+1; k=nz+1
            f1(l, i, j, k)=f1(l, i-cxs(l), j-cys(l), k-czs(l))
         endif
      enddo
   enddo

end subroutine
end module
