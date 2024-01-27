module m_bndapply
contains
subroutine bndapply(f,fbnd,blanking)
   use mod_dimensions
   implicit none
   real, intent(inout) :: f(nx,ny,nz,nl)
   real, intent(in)    :: fbnd(nx,ny,nz,nl)
   logical, intent(in) :: blanking(nx,ny,nz)
   integer i,j,k

!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(f,fbnd)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         if (blanking(i,j,k)) f(i,j,k,1:nl)=fbnd(i,j,k,1:nl)
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO
end subroutine
end module
