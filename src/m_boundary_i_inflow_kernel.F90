module m_boundary_i_inflow_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
subroutine boundary_i_inflow_kernel(f,uvel,rho0,udir)
! Inflow outflow boundary conditions in i-direction.
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions
   use mod_D3Q27setup

   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: uvel(nz)
   real, value   :: udir
   real, value   :: rho0
   integer i,j,k,l,ka
   real, parameter :: pi=3.1415927410125732
   real tmp,wl,cxl,cyl,czl

   i =  1
#ifdef _CUDA
   j = threadIdx%y + (blockIdx%y-1)*blockDim%y-1
   k = threadIdx%z + (blockIdx%z-1)*blockDim%z-1
   if (j > ny+1) return
   if (k > nz+1) return
#else
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(j,k,l,tmp) SHARED(f, cxs, cys, czs, uvel, rho0, udir, weights, ka)
      do k=0,nz+1
      do j=0,ny+1
#endif
! Inflow is formula 5.26 from Kruger 2016
         ka = min(max(k,1), nz)
         do l=1,nl
            wl = weights(l)
            cxl = real(cxs(l)); cyl = real(cys(l)); czl = real(czs(l))
            f(l,0,j,k) = f(l,1,j,k) - 2.0*wl*rho0*(cxl*uvel(ka)*cos(udir*pi/180.0) + cyl*uvel(ka)*sin(udir*pi/180.0))/cs2
         end do


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
