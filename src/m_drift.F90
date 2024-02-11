module m_drift
contains
subroutine drift(f,feq)
! f enter routine in feq following collisions
! f is returned in f
   use mod_dimensions
   use mod_D3Q27setup
   implicit none
   real, intent(out) :: f(nx,ny,nz,nl)
   real, intent(in)  :: feq(nx,ny,nz,nl)
   integer l,i,j,k,i1,j1,k1

!$OMP PARALLEL DO PRIVATE(i,j,k,l,i1,j1,k1) SHARED(f, feq, cxs, cys, czs)
      do l = 1, nl
         do k=1,nz
            k1=mod(nz+k-1-czs(l),nz)+1
            do j=1,ny
               j1=mod(ny+j-1-cys(l),ny)+1
               do i=1,nx
                  i1=mod(nx+i-1-cxs(l),nx)+1
                  f(i,j,k,l) = feq(i1,j1,k1,l)
               enddo
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO

end subroutine
end module
