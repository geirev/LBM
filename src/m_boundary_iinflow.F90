module m_boundary_iinflow
contains
#ifdef _CUDA
   attributes(global)&
#endif
subroutine boundary_iinflow(f,uvel,rho0,udir)
! Inflow outflow boundary conditions in i-direction.
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions
   use mod_D3Q27setup

   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: uvel(nz)
   real, value, intent(in)   :: udir
   real, value, intent(in)   :: rho0
   integer i,j,k,l,ka
   real, parameter :: pi=3.1415927410125732
   real tmp

#ifdef _CUDA
   i = 1
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z
   if (j > nz+2) return
   if (k > nz+2) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(j,k,l,tmp) SHARED(f, nx, ny, nz, nl, cxs, cys, czs, uvel, rho0, udir, weights, ka, pi)
      do k=0,nz+1
      do j=0,ny+1
#endif
! Inflow is formula 5.26 from Kruger 2016
         ka = min(max(k,1), nz)
         do l=1,nl
            f(l,0,j,k) = f(l,1,j,k) - 2.0*weights(l)*rho0*(real(cxs(l))*uvel(ka)*cos(udir*pi/180)&
                                                         + real(cys(l))*uvel(ka)*sin(udir*pi/180)&
                                                         + real(czs(l))*0)/cs2
         enddo

! bounce back specifically coded for current ordering. It could be using a general formula with bounce(l) array
         do l=2,nl-1,2
            tmp=f(l,0,j,k)
            if (cxs(l)==1)   f(l  ,0,j,k)=f(l+1,0,j,k)
            if (cxs(l+1)==1) f(l+1,0,j,k)=tmp
         enddo

! Outflow
         do l=1,nl
            f(l,nx+1,j,k)=f(l,nx,j,k)                                 ! 0th order extrapolation
         enddo
#ifndef _CUDA
      enddo
      enddo
#endif

end subroutine
end module
