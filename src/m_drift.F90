module m_drift
contains
subroutine drift(f,fwork,cxs,cys)
   use mod_dimensions
   implicit none
   real, intent(inout) :: f(nx,ny,nl)
   real, intent(inout) :: fwork(nx,ny,nl)
   integer, intent(in) :: cxs(nl)
   integer, intent(in) :: cys(nl)
   integer l,i,j,i1,j1

!$OMP PARALLEL DO PRIVATE(j,i,l,i1,j1) SHARED(f, fwork, cxs, cys)
      do l = 1, nl
         do j=1,ny
            j1=mod(ny+j-1-cys(l),ny)+1
            do i=1,nx
               fwork(i,j,l) = f(i,j1,l)
            enddo
         enddo
         do j=1,ny
            do i=1,nx
               i1=mod(nx+i-1-cxs(l),nx)+1
               f(i,j,l) = fwork(i1,j,l)
            enddo
         enddo
!         f(:,:,l) = cshift(f(:,:,l),-cxs(l), 1)
!         f(:,:,l) = cshift(f(:,:,l),-cys(l), 2)
      enddo
!$OMP END PARALLEL DO

end subroutine
end module
