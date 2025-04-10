module m_bndbounceback
contains
subroutine bndbounceback(f,blanking)
   use mod_dimensions
!   use m_readinfile,   only : kbndbb
   use m_wtime
   implicit none
   logical, intent(in)    :: blanking(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real tmp
   integer i,j,k,l
   integer, parameter :: icpu=9
   call cpustart()

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, tmp) SHARED(f, blanking)
   do k=0,nz+1
   do j=0,ny+1
   do i=0,nx+1
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
