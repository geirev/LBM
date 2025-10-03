module m_fillghosts
contains

subroutine fillghosts(f)
   use mod_dimensions, only : nx, ny, nz
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(out)     :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   integer i,j,k,i1,j1,k1
#ifdef _CUDA
   attributes(device) :: f
#endif

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, i1, j1, k1) SHARED(f)
#endif
   do k=0,nz+1
   do j=0,ny+1
   do i=0,nx+1
      i1=min(max(i,1),nx)
      j1=min(max(j,1),ny)
      k1=min(max(k,1),nz)
      f(:,i,j,k)=f(:,i1,j1,k1)
   enddo
   enddo
   enddo
end subroutine
end module
