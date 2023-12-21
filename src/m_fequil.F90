module m_fequil
contains
subroutine fequil(feq,rho,u,v,weights,cxs,cys)
   use mod_dimensions
   implicit none
   real,    intent(in) :: rho(nx,ny)
   real,    intent(in) :: u(nx,ny)
   real,    intent(in) :: v(nx,ny)
   real,    intent(in) :: weights(nl)
   integer, intent(in) :: cxs(nl)
   integer, intent(in) :: cys(nl)
   real,    intent(inout):: feq(nx,ny,nl)
   integer i,j,l

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, l) SHARED(feq, rho, u, v, weights, cxs, cys)
   do l=1,nl
   do j=1,ny
   do i=1,nx
      feq(i,j,l) = rho(i,j)*weights(l) *(1.0                                                       &
                                        +3.0*(real(cxs(l))*u(i,j) + real(cys(l))*v(i,j))           &
                                        +9.0*(real(cxs(l))*u(i,j) + real(cys(l))*v(i,j))**2/2.0    &
                                        -3.0*((u(i,j))**2 + (v(i,j))**2)/2.0)
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

end subroutine
end module
