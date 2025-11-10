module m_macrovars_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine macrovars_kernel(f, rho, u, v, w)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl,cxs,cys,czs
   implicit none
   real, intent(in)  :: f(nl, 0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(out) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(out) ::   u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(out) ::   v(0:nx+1,0:ny+1,0:nz+1)
   real, intent(out) ::   w(0:nx+1,0:ny+1,0:nz+1)
   integer :: i, j, k, l
   real tmpf,tmpr,tmpu,tmpv,tmpw,invrho

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x-1
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y-1
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z-1
   if (i < 0 .or. i > nx+1) return
   if (j < 0 .or. j > ny+1) return
   if (k < 0 .or. k > nz+1) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l,tmpf,tmpr,tmpu,tmpv,tmpw,invrho) SHARED(u, v, w, rho, f, cxs, cys, czs)
   do k=0,nz+1
   do j=0,ny+1
   do i=0,nx+1
#endif

      tmpr=0.0; tmpu=0.0; tmpv=0.0; tmpw=0.0

      do l = 1, nl
         tmpf=f(l,i,j,k)
         tmpr=tmpr+tmpf
         tmpu=tmpu + cxs(l)*tmpf
         tmpv=tmpv + cys(l)*tmpf
         tmpw=tmpw + czs(l)*tmpf
      enddo
      invrho=1.0/tmpr

      rho(i,j,k)=tmpr
      u(i,j,k)=tmpu * invrho
      v(i,j,k)=tmpv * invrho
      w(i,j,k)=tmpw * invrho

#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif


end subroutine
end module
