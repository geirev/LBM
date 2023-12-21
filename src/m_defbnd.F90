module m_defbnd
contains
subroutine defbnd(fbnd,f,blanking)
      use mod_dimensions
      implicit none
      logical, intent(in)    :: blanking(nx,ny)
      real,    intent(in)    :: f(nx,ny,nl)
      real,    intent(inout) :: fbnd(nx,ny,nl)
      integer i,j

!!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j) SHARED(fbnd, f, blanking)
      do j=1,ny
      do i=1,nx
         if (blanking(i,j)) then
            fbnd(i,j,1)=f(i,j,1)
            fbnd(i,j,2)=f(i,j,6)
            fbnd(i,j,3)=f(i,j,7)
            fbnd(i,j,4)=f(i,j,8)
            fbnd(i,j,5)=f(i,j,9)
            fbnd(i,j,6)=f(i,j,2)
            fbnd(i,j,7)=f(i,j,3)
            fbnd(i,j,8)=f(i,j,4)
            fbnd(i,j,9)=f(i,j,5)
         endif
      enddo
      enddo
!!$OMP END PARALLEL DO

end subroutine
end module
