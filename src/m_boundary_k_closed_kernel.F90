module m_boundary_k_closed_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_k_closed_kernel(f1,f2,kplane,opt)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl,cxs,cys,czs
   implicit none
   real, intent(inout) :: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   integer, value      :: kplane
   integer, value      :: opt
   integer :: i, j, k, l, m, kghost
#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   if (i < 1 .or. i > nx) return
   if (j < 1 .or. j > ny) return
#endif
   k=1
   if (kplane == 1 ) then
      kghost=0
   elseif (kplane == nz)  then
      kghost=nz+1
   endif

#ifndef _CUDA
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, k, l, m) SHARED(f1,f2, cxs, cys, czs, kghost, kplane, opt)
   do j=1,ny
   do i=1,nx
#endif
      do l = 1, nl
         if ((kplane == 1 .and. czs(l) < 0) .or. (kplane == nz .and. czs(l) > 0)) then
            do m = 1, nl
               if (cxs(m) == opt*cxs(l) .and. cys(m) == opt*cys(l) .and. czs(m) == -czs(l)) then
                 f1(m,i,j,kghost) = f2(l,i,j, kghost)
                 f1(l,i,j,kghost) = f1(l,i-cxs(l),j-cys(l),kplane)
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
