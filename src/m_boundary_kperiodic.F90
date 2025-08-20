module m_boundary_kperiodic
contains
subroutine boundary_kperiodic(f)
   use mod_dimensions
   use mod_D3Q27setup, only : nl

   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
#endif
   integer i,j,l

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do j=0,ny+1
      do i=0,nx+1
      do l=1,nl
         f(l,i,j,0)   =f(l,i,j,nz)
         f(l,i,j,nz+1)=f(l,i,j,1)
      enddo
      enddo
      enddo

end subroutine
end module
