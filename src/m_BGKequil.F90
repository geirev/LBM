module m_BGKequil
contains
subroutine BGKequil(feq,f,rho,u,v,w)
   use mod_dimensions
   use mod_D3Q27setup
   use m_wtime
   implicit none
   real,    intent(in)    :: rho(nx,ny,nz)
   real,    intent(in)    :: u(nx,ny,nz)
   real,    intent(in)    :: v(nx,ny,nz)
   real,    intent(in)    :: w(nx,ny,nz)
   real,    intent(out)   :: feq(0:nx+1,0:ny+1,0:nz+1,nl)
   real,    intent(inout) :: f(0:nx+1,0:ny+1,0:nz+1,nl)
   integer i,j,k,l
   integer, parameter :: icpu=5
   call cpustart()


!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l) SHARED(feq, f, rho, u, v, w, weights, cxs, cys, czs)
   do l=1,nl
      do k=1,nz
      do j=1,ny
      do i=1,nx

! equilibrium distribution
         feq(i,j,k,l) = rho(i,j,k)*weights(l)*(1.0                                                          &
                      +3.0*(real(cxs(l))*u(i,j,k) + real(cys(l))*v(i,j,k) + real(czs(l))*w(i,j,k))          &
                      +9.0*(real(cxs(l))*u(i,j,k) + real(cys(l))*v(i,j,k) + real(czs(l))*w(i,j,k))**2/2.0   &
                      -3.0*((u(i,j,k))**2 + (v(i,j,k))**2 + (w(i,j,k))**2)/2.0)

! non-equilibrium distribution
         f(i,j,k,l)=f(i,j,k,l)-feq(i,j,k,l)

      enddo
      enddo
      enddo
   enddo
!$OMP END PARALLEL DO
   call cpufinish(icpu)

end subroutine
end module
