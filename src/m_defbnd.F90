module m_defbnd
contains
subroutine defbnd(fbnd,f,blanking)
      use mod_dimensions
      implicit none
      logical, intent(in)    :: blanking(nx,ny,nz)
      real,    intent(in)    :: f(nx,ny,nz,nl)
      real,    intent(inout) :: fbnd(nx,ny,nz,nl)
      integer i,j,k,l

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k) SHARED(fbnd, f, blanking)
      do k=1,nz
      do j=1,ny
      do i=1,nx
         if (blanking(i,j,k)) then
            fbnd(i,j,k,1)=f(i,j,k,1)
            do l=2,nl-1,2
               fbnd(i,j,k,l)=f(i,j,k,l+1)
               fbnd(i,j,k,l+1)=f(i,j,k,l)
            enddo
         endif
      enddo
      enddo
      enddo
!$OMP END PARALLEL DO

end subroutine
end module
