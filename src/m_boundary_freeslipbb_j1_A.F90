module m_boundary_freeslipbb_j1_A
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine boundary_freeslipbb_j1_A(f,nx2,ny2,nz2,nl,cxs,cys,czs)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx2, ny2, nz2, nl
   integer             :: cxs(nl),cys(nl),czs(nl)
   real, intent(inout) :: f(nl,nx2,ny2,nz2)
   integer :: i, k, l, m
#ifdef _CUDA
   integer j
   attributes(device) :: f
   attributes(device) :: cxs,cys,czs
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = 1
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 2 .or. i > nx2-1) return
   if (k < 2 .or. k > nz2-1) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, k, l, m) SHARED(f, cxs, cys, czs, nx2, ny2, nz2, nl)
   do k=2,nz2-1
   do i=2,nx2-1
#endif
      do l = 1, nl
         if (cys(l) < 0) then
            do m = 1, nl
               if (cxs(m) == cxs(l) .and. cys(m) == -cys(l) .and. czs(m) ==  czs(l)) then
                 f(m,i,1,k)=f(l,i,1,k)
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
