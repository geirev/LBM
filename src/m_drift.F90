module m_drift
contains
subroutine drift(f,fwork)
   use mod_dimensions
   use mod_D3Q27setup
   implicit none
   real, intent(inout) :: f(nx,ny,nz,nl)
   real, intent(inout) :: fwork(nx,ny,nz,nl)
   integer l,i,j,k,i1,j1,k1

!$OMP PARALLEL DO PRIVATE(i,j,k,l,i1,j1,k1) SHARED(f, fwork, cxs, cys, czs)
      do l = 1, nl
         do k=1,nz
            do j=1,ny
               do i=1,nx
                  i1=mod(nx+i-1-cxs(l),nx)+1
                  fwork(i,j,k,l) = f(i1,j,k,l)
               enddo
            enddo
         enddo

         do k=1,nz
            do j=1,ny
               j1=mod(ny+j-1-cys(l),ny)+1
               f(:,j,k,l) = fwork(:,j1,k,l)
            enddo
         enddo

         do k=1,nz
            k1=mod(nz+k-1-czs(l),nz)+1
            fwork(:,:,k,l) = f(:,:,k1,l)
         enddo

         f(:,:,:,l) = fwork(:,:,:,l)

!         f(:,:,:,l) = cshift(f(:,:,l),-cxs(l), 1)
!         f(:,:,:,l) = cshift(f(:,:,l),-cys(l), 2)
!         f(:,:,:,l) = cshift(f(:,:,l),-czs(l), 3)
      enddo
!$OMP END PARALLEL DO

end subroutine
end module
