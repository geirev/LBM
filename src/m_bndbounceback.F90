module m_bndbounceback
contains
subroutine bndbounceback(f,blanking)
   use mod_dimensions
   use m_wtime
   implicit none
   logical, intent(in)    :: blanking(nx,ny,nz)
   real,    intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real tmp
   integer i,j,k,l
   integer, parameter :: icpu=9
   call cpustart()

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, tmp) SHARED(f, blanking)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      if (blanking(i,j,k)) then
         do l=2,nl-1,2
            tmp=f(l,i,j,k)
            f(l,i,j,k)=f(l+1,i,j,k)
            f(l+1,i,j,k)=tmp
         enddo
      endif
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

   call cpufinish(icpu)
end subroutine
end module
