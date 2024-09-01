module m_applyturbines
! based on the Kupershtokh (2004) method (see Eq 6.31 in Kruger book).
contains
subroutine applyturbines(f,df)
   use mod_dimensions
   use m_readinfile
   use m_wtime
   implicit none
   real, intent(inout) :: f(0:nx+1,0:ny+1,0:nz+1,nl) ! distribution
   real, intent(in)    :: df(-ieps:ieps,ny,nz,nl,nturbines)     ! forcing distributions
   integer n,ip
   integer, parameter :: icpu=7
   call cpustart()
   do n=1,nturbines
      ip=ipos(n)
      f(ip-ieps:ip+ieps,1:ny,1:nz,1:nl) = f(ip-ieps:ip+ieps,1:ny,1:nz,1:nl) + df(-ieps:ieps,1:ny,1:nz,1:nl,n)
   enddo
   call cpufinish(icpu)
end subroutine
end module
