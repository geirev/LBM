module m_macrovars_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine macrovars_kernel(f, rho, u, v, w, nx, ny, nz, nl)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value :: nx, ny, nz, nl
   real, intent(in)  :: f(nl, 0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(out) :: rho(nx,ny,nz)
   real, intent(out) :: u(nx,ny,nz)
   real, intent(out) :: v(nx,ny,nz)
   real, intent(out) :: w(nx,ny,nz)
   integer :: i, j, k, l
   real fl,tmpr,tmpu,tmpv,tmpw,inv_rho

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx) return
   if (j > ny) return
   if (k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(u, v, w, rho, f)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
       rho(i,j,k)=0.0
       do l = 1, nl
          rho(i,j,k)=rho(i,j,k)+f(l,i,j,k)
       enddo

       u(i,j,k) = f( 2,i,j,k)  - f( 3,i,j,k) + f( 8,i,j,k) - f( 9,i,j,k) + f(10,i,j,k) - f(11,i,j,k) - f(12,i,j,k) &
                + f(13,i,j,k)  - f(16,i,j,k) + f(17,i,j,k)
#ifndef D3Q19
       u(i,j,k) = u(i,j,k)     - f(20,i,j,k) + f(21,i,j,k) - f(22,i,j,k) + f(23,i,j,k) + f(24,i,j,k) - f(25,i,j,k) &
                               - f(26,i,j,k) + f(27,i,j,k)
#endif
       u(i,j,k) = u(i,j,k)/rho(i,j,k)

       v(i,j,k) = f( 4,i,j,k)  - f( 5,i,j,k) + f( 8,i,j,k) - f( 9,i,j,k) - f(10,i,j,k) + f(11,i,j,k) + f(14,i,j,k) &
                               - f(15,i,j,k) - f(18,i,j,k) + f(19,i,j,k)
#ifndef D3Q19
       v(i,j,k) = v(i,j,k)     + f(20,i,j,k) - f(21,i,j,k) - f(22,i,j,k) + f(23,i,j,k) + f(24,i,j,k) - f(25,i,j,k) &
                               + f(26,i,j,k) - f(27,i,j,k)
#endif
       v(i,j,k) = v(i,j,k)/rho(i,j,k)

       w(i,j,k) = -f( 6,i,j,k) + f( 7,i,j,k) - f(12,i,j,k) + f(13,i,j,k) + f(14,i,j,k) - f(15,i,j,k) + f(16,i,j,k) &
                               - f(17,i,j,k) + f(18,i,j,k) - f(19,i,j,k)
#ifndef D3Q19
       w(i,j,k) = w(i,j,k)     + f(20,i,j,k) - f(21,i,j,k) - f(22,i,j,k) + f(23,i,j,k) - f(24,i,j,k) + f(25,i,j,k) &
                               - f(26,i,j,k) + f(27,i,j,k)
#endif
       w(i,j,k) =  w(i,j,k)/rho(i,j,k)
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif


end subroutine
end module




!!  module m_macrovars_kernel
!!  contains
!!  #ifdef _CUDA
!!     attributes(global)&
!!  #endif
!!     subroutine macrovars_kernel(f, rho, u, v, w, nx, ny, nz, nl)
!!  #ifdef _CUDA
!!     use cudafor
!!  #endif
!!     implicit none
!!     integer, value :: nx, ny, nz, nl
!!     real, intent(in)  :: f(nl, 0:nx+1, 0:ny+1, 0:nz+1)
!!     real, intent(out) :: rho(nx,ny,nz)
!!     real, intent(out) :: u(nx,ny,nz)
!!     real, intent(out) :: v(nx,ny,nz)
!!     real, intent(out) :: w(nx,ny,nz)
!!  
!!     integer :: i, j, k
!!     ! locals for density and velocities
!!     real :: tmp_rho, tmp_u, tmp_v, tmp_w, inv_rho
!!     ! temporaries for all directions
!!     real :: f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13
!!     real :: f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26,f27
!!  
!!  #ifdef _CUDA
!!     i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
!!     j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
!!     k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
!!     if (i < 1 .or. i > nx) return
!!     if (j < 1 .or. j > ny) return
!!     if (k < 1 .or. k > nz) return
!!  #else
!!  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k, &
!!  !$OMP    f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13, &
!!  !$OMP    f14,f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26,f27, &
!!  !$OMP    tmp_rho,tmp_u,tmp_v,tmp_w,inv_rho) &
!!  !$OMP SHARED(f,rho,u,v,w,nx,ny,nz)
!!     do k=1,nz
!!     do j=1,ny
!!     do i=1,nx
!!  #endif
!!  
!!        ! --- Load once into registers ---
!!        f1  = f( 1,i,j,k);  f2  = f( 2,i,j,k);  f3  = f( 3,i,j,k)
!!        f4  = f( 4,i,j,k);  f5  = f( 5,i,j,k);  f6  = f( 6,i,j,k)
!!        f7  = f( 7,i,j,k);  f8  = f( 8,i,j,k);  f9  = f( 9,i,j,k)
!!        f10 = f(10,i,j,k);  f11 = f(11,i,j,k);  f12 = f(12,i,j,k)
!!        f13 = f(13,i,j,k);  f14 = f(14,i,j,k);  f15 = f(15,i,j,k)
!!        f16 = f(16,i,j,k);  f17 = f(17,i,j,k);  f18 = f(18,i,j,k)
!!        f19 = f(19,i,j,k);  f20 = f(20,i,j,k);  f21 = f(21,i,j,k)
!!        f22 = f(22,i,j,k);  f23 = f(23,i,j,k);  f24 = f(24,i,j,k)
!!        f25 = f(25,i,j,k);  f26 = f(26,i,j,k);  f27 = f(27,i,j,k)
!!  
!!        ! --- Density ---
!!        tmp_rho = f1+f2+f3+f4+f5+f6+f7+f8+f9+f10+f11+f12+f13+ &
!!                  f14+f15+f16+f17+f18+f19+f20+f21+f22+f23+f24+f25+f26+f27
!!        rho(i,j,k) = tmp_rho
!!  
!!        ! --- Velocities ---
!!        tmp_u =  f2 - f3 + f8 - f9 + f10 - f11 - f12 + f13 - f16 + f17
!!  #ifndef D3Q19
!!        tmp_u = tmp_u - f20 + f21 - f22 + f23 + f24 - f25 - f26 + f27
!!  #endif
!!  
!!        tmp_v =  f4 - f5 + f8 - f9 - f10 + f11 + f14 - f15 - f18 + f19
!!  #ifndef D3Q19
!!        tmp_v = tmp_v + f20 - f21 - f22 + f23 + f24 - f25 + f26 - f27
!!  #endif
!!  
!!        tmp_w = -f6 + f7 - f12 + f13 + f14 - f15 + f16 - f17 + f18 - f19
!!  #ifndef D3Q19
!!        tmp_w = tmp_w + f20 - f21 - f22 + f23 - f24 + f25 - f26 + f27
!!  #endif
!!  
!!        inv_rho = 1.0 / tmp_rho
!!        u(i,j,k) = tmp_u * inv_rho
!!        v(i,j,k) = tmp_v * inv_rho
!!        w(i,j,k) = tmp_w * inv_rho
!!  
!!  #ifndef _CUDA
!!     enddo
!!     enddo
!!     enddo
!!  !$OMP END PARALLEL DO
!!  #endif
!!  
!!     end subroutine
!!  end module

!!  module m_macrovars_kernel
!!  contains
!!  #ifdef _CUDA
!!     attributes(global)&
!!  #endif
!!     subroutine macrovars_kernel(f, rho, u, v, w, cxs, cys, czs, nx2, ny2, nz2, nl)
!!  #ifdef _CUDA
!!     use cudafor
!!  #endif
!!     implicit none
!!     integer, value :: nx2, ny2, nz2, nl
!!     real, intent(in)  :: f(nl, nx2, ny2, nz2)
!!     integer, intent(in)  :: cxs(nl)
!!     integer, intent(in)  :: cys(nl)
!!     integer, intent(in)  :: czs(nl)
!!     real, intent(out) :: rho(nx2-2,ny2-2,nz2-2)
!!     real, intent(out) :: u(nx2-2,ny2-2,nz2-2)
!!     real, intent(out) :: v(nx2-2,ny2-2,nz2-2)
!!     real, intent(out) :: w(nx2-2,ny2-2,nz2-2)
!!     integer :: i, j, k, l
!!  #ifdef _CUDA
!!     attributes(device) :: rho
!!     attributes(device) :: u
!!     attributes(device) :: v
!!     attributes(device) :: w
!!     attributes(device) :: f
!!     attributes(device) :: cxs
!!     attributes(device) :: cys
!!     attributes(device) :: czs
!!  #endif
!!  
!!  
!!  #ifdef _CUDA
!!     i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
!!     j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
!!     k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
!!     if (i > nx2-2) return
!!     if (j > ny2-2) return
!!     if (k > nz2-2) return
!!  #else
!!  !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(u, v, w, rho, f, cxs, cys, czs, nz2, ny2, nx2, nl )
!!     do k=1,nz2-2
!!     do j=1,ny2-2
!!     do i=1,nx2-2
!!  #endif
!!  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        rho(i,j,k)=f(1,i+1,j+1,k+1)
!!        u(i,j,k)=0.0
!!        v(i,j,k)=0.0
!!        w(i,j,k)=0.0
!!        do l = 2, nl
!!           rho(i,j,k)=rho(i,j,k)+             f(l,i+1,j+1,k+1)
!!             u(i,j,k)=  u(i,j,k)+real(cxs(l))*f(l,i+1,j+1,k+1)
!!             v(i,j,k)=  v(i,j,k)+real(cys(l))*f(l,i+1,j+1,k+1)
!!             w(i,j,k)=  w(i,j,k)+real(czs(l))*f(l,i+1,j+1,k+1)
!!        enddo
!!        u(i,j,k)=u(i,j,k)/rho(i,j,k)
!!        v(i,j,k)=v(i,j,k)/rho(i,j,k)
!!        w(i,j,k)=w(i,j,k)/rho(i,j,k)
!!  #ifndef _CUDA
!!      enddo
!!      enddo
!!      enddo
!!  !$OMP END PARALLEL DO
!!  #endif
!!  
!!  
!!  
!!  end subroutine
!!  end module
