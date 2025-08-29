module m_boundary_iperiodic
contains
subroutine boundary_iperiodic(f)
   use mod_dimensions
   use mod_D3Q27setup, only : nl

   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
#endif

   integer j,k,l
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
      do l=1,nl
         f(l,0,j,k)   =f(l,nx,j,k)
         f(l,nx+1,j,k)=f(l,1,j,k)
      enddo
      enddo
      enddo

end subroutine
end module
