module m_bndbounceback
contains
subroutine bndbounceback(f,blanking)
   use mod_dimensions
   implicit none
   logical, intent(in)    :: blanking(nx,ny,nz)
   real,    intent(inout) :: f(0:nx+1,ny,nz,nl)
   real tmp
   integer i,j,k,l

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, tmp) SHARED(f, blanking)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (blanking(i,j,k)) then
         do l=2,nl-1,2
            tmp=f(i,j,k,l)
            f(i,j,k,l)=f(i,j,k,l+1)
            f(i,j,k,l+1)=tmp
         enddo
      endif
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

end subroutine
end module
