module m_boundary_closed_corners
contains
subroutine boundary_closed_corners(f1,f2,opt)
   use mod_dimensions,    only: nx,ny,nz
   use mod_D3Q27setup,    only: cxs,cys,czs,nl
   implicit none
   real, intent(inout) :: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, value      :: opt   ! -1 = noslip, +1 = freeslip

#ifdef _CUDA
   attributes(device) :: f1
   attributes(device) :: f2
#endif

   integer :: i,j,k,l,m
   integer :: ic,jc,kc      ! corner coordinates
   integer :: iin,jin,kin   ! interior neighbor

   !===========================================================
   ! 8 CORNERS:
   !   (ic,jc,kc) ∈ {0,nx+1} × {0,ny+1} × {0,nz+1}
   !===========================================================
   integer, dimension(8,3) :: corner = reshape([ &
         0,     0,     0,     & ! corner 1
         0,     0,     nz+1,  & ! corner 2
         0,     ny+1,  0,     & ! corner 3
         0,     ny+1,  nz+1,  & ! corner 4
         nx+1,  0,     0,     & ! corner 5
         nx+1,  0,     nz+1,  & ! corner 6
         nx+1,  ny+1,  0,     & ! corner 7
         nx+1,  ny+1,  nz+1   & ! corner 8
      ], [8,3])

   !-------------------------------
   ! Loop through all 8 corners
   !-------------------------------
   do m=1,8

      ic = corner(m,1)
      jc = corner(m,2)
      kc = corner(m,3)

      ! =======================================================
      ! PASS 1: Set incoming populations at the corner (bounce)
      ! =======================================================

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do l=1,nl

         ! Directions pointing OUT of the fluid domain:
         ! sign(cxs), sign(cys), sign(czs) must point *towards* the corner

         if ( ( (ic==0    .and. cxs(l)<0)  .or. (ic==nx+1 .and. cxs(l)>0) ) .and. &
              ( (jc==0    .and. cys(l)<0)  .or. (jc==ny+1 .and. cys(l)>0) ) .and. &
              ( (kc==0    .and. czs(l)<0)  .or. (kc==nz+1 .and. czs(l)>0) ) ) then

            ! Find mirrored direction m
            do i=1,nl
               if (cxs(i) == opt*cxs(l) .and. &
                   cys(i) == opt*cys(l) .and. &
                   czs(i) == opt*czs(l)) then

                  f1(i,ic,jc,kc) = f2(l,ic,jc,kc)
                  exit
               endif
            enddo

         endif
      enddo


      ! =======================================================
      ! PASS 2: outgoing directions: copy from interior neighbor
      ! =======================================================

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do l=1,nl

         ! Same outgoing criteria as pass 1
         if ( ( (ic==0    .and. cxs(l)<0)  .or. (ic==nx+1 .and. cxs(l)>0) ) .and. &
              ( (jc==0    .and. cys(l)<0)  .or. (jc==ny+1 .and. cys(l)>0) ) .and. &
              ( (kc==0    .and. czs(l)<0)  .or. (kc==nz+1 .and. czs(l)>0) ) ) then

            ! interior neighbor index
            iin = ic - cxs(l)
            jin = jc - cys(l)
            kin = kc - czs(l)

            f1(l,ic,jc,kc) = f1(l, iin, jin, kin)
         endif

      enddo

   enddo

end subroutine
end module

