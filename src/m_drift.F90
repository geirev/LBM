module m_drift
contains
subroutine drift(f,feq)
! f enter routine in feq following collisions
! f is returned in f
   use mod_dimensions
   use mod_D3Q27setup
   use m_wtime
   implicit none
   real, intent(out) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)  :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   integer l,i,j,k,i1,j1,k1
   integer, parameter :: icpu=10
   call cpustart()


!$OMP PARALLEL DO PRIVATE(i,j,k,l,i1,j1,k1) SHARED(f, feq, cxs, cys, czs)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      do l=1,nl
         f(l,i,j,k) = feq(l,i-cxs(l),j-cys(l),k-czs(l))
      enddo
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

    call cpufinish(icpu)

end subroutine
end module
