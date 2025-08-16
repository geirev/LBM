module m_macrovars_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine macrovars_kernel(f, rho, u, v, w, cxs, cys, czs, nx2, ny2, nz2, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value :: nx2, ny2, nz2, nl
   real, intent(in)  :: f(nl, nx2, ny2, nz2)
   integer, intent(in)  :: cxs(nl)
   integer, intent(in)  :: cys(nl)
   integer, intent(in)  :: czs(nl)
   real, intent(out) :: rho(nx2-2,ny2-2,nz2-2)
   real, intent(out) :: u(nx2-2,ny2-2,nz2-2)
   real, intent(out) :: v(nx2-2,ny2-2,nz2-2)
   real, intent(out) :: w(nx2-2,ny2-2,nz2-2)
   integer :: i, j, k, l
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: f
   attributes(device) :: cxs
   attributes(device) :: cys
   attributes(device) :: czs
#endif


#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx2-2) return
   if (j > ny2-2) return
   if (k > nz2-2) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(u, v, w, rho, f, cxs, cys, czs, nz2, ny2, nx2, nl )
   do k=1,nz2-2
   do j=1,ny2-2
   do i=1,nx2-2
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rho(i,j,k)=f(1,i+1,j+1,k+1)
      u(i,j,k)=0.0
      v(i,j,k)=0.0
      w(i,j,k)=0.0
      do l = 2, nl
         rho(i,j,k)=rho(i,j,k)+             f(l,i+1,j+1,k+1)
           u(i,j,k)=  u(i,j,k)+real(cxs(l))*f(l,i+1,j+1,k+1)
           v(i,j,k)=  v(i,j,k)+real(cys(l))*f(l,i+1,j+1,k+1)
           w(i,j,k)=  w(i,j,k)+real(czs(l))*f(l,i+1,j+1,k+1)
      enddo
      u(i,j,k)=u(i,j,k)/rho(i,j,k)
      v(i,j,k)=v(i,j,k)/rho(i,j,k)
      w(i,j,k)=w(i,j,k)/rho(i,j,k)
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif



end subroutine
end module
