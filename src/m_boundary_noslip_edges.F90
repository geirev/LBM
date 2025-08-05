module m_boundary_noslip_edges
contains
subroutine boundary_noslip_edges(f)
   use mod_dimensions
   use mod_D3Q27setup, only : cxs,cys,czs
   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
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
               if (cxs(m) == -cxs(l) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
                  f(m,i,0,0)=f(l,i,0,0)
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
            f(l, i, j, k)=f(l, i-cxs(l), j-cys(l), k-czs(l))
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
               if (cxs(m) == -cxs(l) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
                  f(m,i,ny+1,0)=f(l,i,ny+1,0)
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
            f(l, i, j, k)=f(l, i-cxs(l), j-cys(l), k-czs(l))
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
               if (cxs(m) == -cxs(l) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
                  f(m,i,0,nz+1)=f(l,i,0,nz+1)
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
            f(l, i, j, k)=f(l, i-cxs(l), j-cys(l), k-czs(l))
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
               if (cxs(m) == -cxs(l) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
                  f(m,i,ny+1,nz+1)=f(l,i,ny+1,nz+1)
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
            f(l, i, j, k)=f(l, i-cxs(l), j-cys(l), k-czs(l))
         endif
      enddo
   enddo

end subroutine
end module
