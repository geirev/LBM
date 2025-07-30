module m_boundary_jperiodic
contains
subroutine boundary_jperiodic(f)
   use mod_dimensions

   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
#endif
   integer i,k,l

#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do i=0,nx+1
      do l=1,nl
         f(l,i,0,k)   =f(l,i,ny,k)
         f(l,i,ny+1,k)=f(l,i,1,k)
      enddo
      enddo
      enddo

end subroutine
end module
