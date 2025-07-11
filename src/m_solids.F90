module m_solids
contains
subroutine solids(f,blanking)
   use mod_dimensions
!   use m_readinfile,   only : kbndbb
   use m_wtime
   implicit none
   logical, intent(in)    :: blanking(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: blanking
#endif
   real tmp
   integer i,j,k,l
   integer, parameter :: icpu=9
   call cpustart()

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, tmp) SHARED(f, blanking)
#endif
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
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif

   call cpufinish(icpu)
end subroutine
end module
