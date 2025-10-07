module m_boundary_i_inflow_edges
contains
subroutine boundary_i_inflow_edges(f1,f2)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(inout):: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout):: f2(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f1
   attributes(device) :: f2
#endif
   integer j,k

!inflow
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
         f1(:,0,0,k)=f1(:,1,0,k)
         f1(:,0,ny+1,k)=f1(:,1,ny+1,k)
      enddo
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do j=0,ny+1
         f1(:,0,j,0)=f1(:,1,j,0)
         f1(:,0,j,nz+1)=f1(:,1,j,nz+1)
      enddo

!outflow
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1,ny+1
          f1(:,nx+1,j,k)=f1(:,nx,j,k)
      enddo
      enddo
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1,nz+1
      do j=0,ny+1
          f1(:,nx+1,j,k)=f1(:,nx,j,k)
      enddo
      enddo

end subroutine
end module
