module m_fequil
contains
subroutine fequil(feq,rho,u,v,w)
   use mod_dimensions
   use mod_D3Q27setup
   implicit none
   real,    intent(in) :: rho(nx,ny,nz)
   real,    intent(in) :: u(nx,ny,nz)
   real,    intent(in) :: v(nx,ny,nz)
   real,    intent(in) :: w(nx,ny,nz)
   real,    intent(inout):: feq(0:nx+1,ny,nz,nl)
   integer i,j,k,l

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l) SHARED(feq, rho, u, v, w, weights, cxs, cys, czs)
   do l=1,nl
   do k=1,nz
   do j=1,ny
   do i=1,nx
      feq(i,j,k,l) = rho(i,j,k)*weights(l) *(1.0                                                  &
                     +3.0*(real(cxs(l))*u(i,j,k) + real(cys(l))*v(i,j,k) + real(czs(l))*w(i,j,k)) &
                     +9.0*(real(cxs(l))*u(i,j,k) + real(cys(l))*v(i,j,k) + real(czs(l))*w(i,j,k))**2/2.0    &
                     -3.0*((u(i,j,k))**2 + (v(i,j,k))**2 + (w(i,j,k))**2)/2.0)
   enddo
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

end subroutine
end module
