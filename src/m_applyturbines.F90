module m_applyturbines
contains
subroutine applyturbines(f,df)
   use mod_dimensions
   use m_readinfile
   implicit none
   real, intent(inout) :: f(0:nx+1,0:ny+1,0:nz+1,nl) ! distribution
   real, intent(in)    :: df(ny,nz,nl,nturbines)     ! forcing distributions
   integer n
   do n=1,nturbines
      f(ipos(n),1:ny,1:nz,1:nl) = f(ipos(n),1:ny,1:nz,1:nl) + df(1:ny,1:nz,1:nl,n)
   enddo
end subroutine
end module
