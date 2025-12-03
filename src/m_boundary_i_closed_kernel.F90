module m_boundary_i_closed_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine boundary_i_closed_kernel(f1,f2,iplane,opt)
#ifdef _CUDA
      use cudafor
#endif
      use mod_dimensions, only : nx,ny,nz
      use mod_D3Q27setup, only : nl, cxs, cys, czs
      implicit none

      real, intent(inout) :: f1(nl,0:nx+1,0:ny+1,0:nz+1)
      real, intent(inout) :: f2(nl,0:nx+1,0:ny+1,0:nz+1)
      integer, value      :: iplane      ! =1 or nx
      integer, value      :: opt         ! -1 no-slip, +1 free-slip

      integer :: i, j, k, l, m, ighost

#ifdef _CUDA
      j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
      k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
      if (j < 1 .or. j > ny) return
      if (k < 1 .or. k > nz) return
#endif

      !----------------------------------------------
      ! Determine ghost index
      !----------------------------------------------
      if (iplane == 1) then
         ighost = 0
      elseif (iplane == nx) then
         ighost = nx+1
      endif

#ifndef _CUDA
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(j,k,l,m) SHARED(f1,f2,cxs,cys,czs,iplane,ighost,opt)
      do k=1,nz
      do j=1,ny
#endif


         do l=1,nl
            ! Select populations pointing OUT of fluid toward ghost
            if ( (iplane == 1  .and. cxs(l) < 0) .or. (iplane == nx .and. cxs(l) > 0) ) then

               ! Find reflected direction m
               do m=1,nl
                  if (cxs(m) == -cxs(l) .and. cys(m) ==  opt*cys(l) .and. czs(m) ==  opt*czs(l)) then
                     f1(m,ighost,j,k) = f2(l,ighost,j,k)
                     f1(l,ighost,j,k) = f1(l,iplane, j-cys(l), k-czs(l))
                     exit
                  endif
               enddo

            endif
         enddo

#ifndef _CUDA
      enddo
      enddo
!$OMP END PARALLEL DO
#endif

   end subroutine
end module

