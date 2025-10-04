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
   real, intent(out) :: rho(nx,ny,nz)
   real, intent(out) :: u(nx,ny,nz)
   real, intent(out) :: v(nx,ny,nz)
   real, intent(out) :: w(nx,ny,nz)
   integer :: i, j, k, l
   real tmpf,tmpr,tmpu,tmpv,tmpw,invrho

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i < 1 .or. i > nx) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(u, v, w, rho, f, nx, ny, nz, nl)
   do k=1,nz
   do j=1,ny
   do i=1,nx
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
