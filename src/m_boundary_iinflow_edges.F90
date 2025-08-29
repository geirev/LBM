module m_boundary_iinflow_edges
contains
subroutine boundary_iinflow_edges(f)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
#endif
   integer j,k

!inflow
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
         f(:,0,0,k)=f(:,1,0,k)
         f(:,0,ny+1,k)=f(:,1,ny+1,k)
      enddo
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do j=0,ny+1
         f(:,0,j,0)=f(:,1,j,0)
         f(:,0,j,nz+1)=f(:,1,j,nz+1)
      enddo

!outflow
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1,ny+1
          f(:,nx+1,j,k)=f(:,nx,j,k)
      enddo
      enddo
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1,nz+1
      do j=0,ny+1
          f(:,nx+1,j,k)=f(:,nx,j,k)
      enddo
      enddo

end subroutine
end module
