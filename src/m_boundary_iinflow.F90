module m_boundary_iinflow
contains
subroutine boundary_iinflow(f,uvel)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile, only : rho0,udir

   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: uvel(nz)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uvel
#endif
   integer j,k,l,ka
   real, parameter :: pi=3.1415927410125732
   real tmp

! Inflow outflow boundary conditions in i-direction.
! Inflow is formula 5.26 from Kruger 2016
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
         ka = min(max(k,1), nz)
         do l=1,nl
            f(l,0,j,k) = f(l,1,j,k) - 2.0*weights(l)*rho0*(real(cxs(l))*uvel(ka)*cos(udir*pi/180)&
                                                         + real(cys(l))*uvel(ka)*sin(udir*pi/180)&
                                                         + real(czs(l))*0)/cs2
         enddo
      enddo
      enddo

! bounce back specifically coded for current ordering. It could be using a general formula with bounce(l) array
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
         do l=2,nl-1,2
            tmp=f(l,0,j,k)
            if (cxs(l)==1)   f(l  ,0,j,k)=f(l+1,0,j,k)
            if (cxs(l+1)==1) f(l+1,0,j,k)=tmp
         enddo
      enddo
      enddo

! Outflow
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
!         f(:,nx+1,j,k) = 2.0*f(:,nx,j,k) - f(:,nx-1,j,k)
          f(:,nx+1,j,k)=f(:,nx,j,k)                                 ! 0th order extrapolation
      enddo
      enddo
!  f(:,nx+1,j,k)=f(:,nx,j,k)                                 ! 0th order extrapolation
!  f(:,i,j,k) = (1 - alpha)*f(:,i,j,k) + alpha*feq(:,i,j,k)  ! non-reflective sponge layer

end subroutine
end module
