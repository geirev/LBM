module m_boundary_j_closed_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_j_closed_kernel(f1,f2,jplane,opt)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl,cxs,cys,czs
   implicit none
   real, intent(inout) :: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, value      :: jplane
   integer, value      :: opt
   integer :: i, j, k, l, m, jghost
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 1 .or. i > nx) return
   if (k < 1 .or. k > nz) return
#endif
   if (jplane == 1 ) then
      jghost=0
   elseif (jplane == ny)  then
      jghost=ny+1
   endif

#ifndef _CUDA
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, k, l, m) SHARED(f1,f2, cxs, cys, czs, jplane, jghost, opt)
   do k=1,nz
   do i=1,nx
#endif
      do l = 1, nl
         if ((jplane == 1 .and. cys(l) < 0) .or. (jplane == ny .and. cys(l) > 0)) then
            ! setting both incoming (m) and outgoing (l) populations at the ghost node.
            do m = 1, nl
               if (cxs(m) == opt*cxs(l) .and. cys(m) == -cys(l) .and. czs(m) ==  opt*czs(l)) then
                 f1(m,i,jghost,k)=f2(l,i,jghost,k)
                 f1(l,i,jghost,k) = f1(l, i-cxs(l), jplane, k-czs(l))
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
