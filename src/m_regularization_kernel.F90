module m_regularization_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine regularization_kernel(f, feq, u, v, w, nx, ny, nz, nl, h2, h3, weights, inv2cs4, inv6cs6, ihrr)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx, ny, nz, nl
   real, intent(inout) :: f(nl,nx+2,ny+2,nz+2)
   real, intent(inout) :: feq(nl,nx+2,ny+2,nz+2)
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)
   real, intent(in)    :: h2(3,3,nl)
   real, intent(in)    :: h3(3,3,3,nl)
   real, intent(in)    :: weights(nl)
   real, value         :: inv2cs4
   real, value         :: inv6cs6
   integer, value      :: ihrr


   real :: vel(3)
   real :: a1_2(3,3)
   real :: a1_3(3,3,3)

   integer :: i, j, k, l, p, q, r, i1, j1, k1
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
   attributes(device) :: h2
   attributes(device) :: h3
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: weights
   i = threadidx%x + (blockidx%x - 1) * blockdim%x
   j = threadidx%y + (blockidx%y - 1) * blockdim%y
   k = threadidx%z + (blockidx%z - 1) * blockdim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(i, j, k, l, p, q, r, i1, j1, k1, vel, a1_2, a1_3)&
!$OMP             & SHARED(f, feq, u, v, w, nx, ny, nz, nl, h2, h3, weights, inv2cs4, inv6cs6, ihrr)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      i1=i+1
      j1=j+1
      k1=k+1

!      do l=1,nl
!          f(l, i1, j1, k1) = f(l, i1, j1, k1) - feq(l, i1, j1, k1)
!      enddo
      f( 1, i1, j1, k1) = f( 1, i1, j1, k1) - feq( 1, i1, j1, k1)
      f( 2, i1, j1, k1) = f( 2, i1, j1, k1) - feq( 2, i1, j1, k1)
      f( 3, i1, j1, k1) = f( 3, i1, j1, k1) - feq( 3, i1, j1, k1)
      f( 4, i1, j1, k1) = f( 4, i1, j1, k1) - feq( 4, i1, j1, k1)
      f( 5, i1, j1, k1) = f( 5, i1, j1, k1) - feq( 5, i1, j1, k1)
      f( 6, i1, j1, k1) = f( 6, i1, j1, k1) - feq( 6, i1, j1, k1)
      f( 7, i1, j1, k1) = f( 7, i1, j1, k1) - feq( 7, i1, j1, k1)
      f( 8, i1, j1, k1) = f( 8, i1, j1, k1) - feq( 8, i1, j1, k1)
      f( 9, i1, j1, k1) = f( 9, i1, j1, k1) - feq( 9, i1, j1, k1)
      f(10, i1, j1, k1) = f(10, i1, j1, k1) - feq(10, i1, j1, k1)
      f(11, i1, j1, k1) = f(11, i1, j1, k1) - feq(11, i1, j1, k1)
      f(12, i1, j1, k1) = f(12, i1, j1, k1) - feq(12, i1, j1, k1)
      f(13, i1, j1, k1) = f(13, i1, j1, k1) - feq(13, i1, j1, k1)
      f(14, i1, j1, k1) = f(14, i1, j1, k1) - feq(14, i1, j1, k1)
      f(15, i1, j1, k1) = f(15, i1, j1, k1) - feq(15, i1, j1, k1)
      f(16, i1, j1, k1) = f(16, i1, j1, k1) - feq(16, i1, j1, k1)
      f(17, i1, j1, k1) = f(17, i1, j1, k1) - feq(17, i1, j1, k1)
      f(18, i1, j1, k1) = f(18, i1, j1, k1) - feq(18, i1, j1, k1)
      f(19, i1, j1, k1) = f(19, i1, j1, k1) - feq(19, i1, j1, k1)
#ifndef D3Q19
      f(20, i1, j1, k1) = f(20, i1, j1, k1) - feq(20, i1, j1, k1)
      f(21, i1, j1, k1) = f(21, i1, j1, k1) - feq(21, i1, j1, k1)
      f(22, i1, j1, k1) = f(22, i1, j1, k1) - feq(22, i1, j1, k1)
      f(23, i1, j1, k1) = f(23, i1, j1, k1) - feq(23, i1, j1, k1)
      f(24, i1, j1, k1) = f(24, i1, j1, k1) - feq(24, i1, j1, k1)
      f(25, i1, j1, k1) = f(25, i1, j1, k1) - feq(25, i1, j1, k1)
      f(26, i1, j1, k1) = f(26, i1, j1, k1) - feq(26, i1, j1, k1)
      f(27, i1, j1, k1) = f(27, i1, j1, k1) - feq(27, i1, j1, k1)
#endif

      if (ihrr == 1) then
! copy u,v,w to vel(1:3)
         vel(1)=u(i,j,k)
         vel(2)=v(i,j,k)
         vel(3)=w(i,j,k)

! computing a1_2
            l=1
            do q=1,3
            do p=1,3
               a1_2(p,q) = h2(p,q,l)*f(l,i1,j1,k1)
            enddo
            enddo

            do l=2,nl
               do q=1,3
               do p=1,3
                  a1_2(p,q) = a1_2(p,q) + h2(p,q,l)*f(l,i1,j1,k1)
               enddo
               enddo
            enddo


! computing a1_3
            do r=1,3
            do q=1,3
            do p=1,3
               a1_3(p,q,r)=vel(p)*a1_2(q,r) + vel(q)*a1_2(r,p) +  vel(r)*a1_2(p,q)
            enddo
            enddo
            enddo


! scale a1_2 and a1_3
            do q=1,3
            do p=1,3
               a1_2(p,q) = a1_2(p,q)*inv2cs4
            enddo
            enddo

            do r=1,3
            do q=1,3
            do p=1,3
               a1_3(p,q,r) = a1_3(p,q,r)*inv6cs6
            enddo
            enddo
            enddo


! second order
         do l=1,nl
            f(l,i1,j1,k1)=0.0
            do q=1,3
            do p=1,3
               f(l,i1,j1,k1)=f(l,i1,j1,k1) + h2(p,q,l)*a1_2(p,q)
            enddo
            enddo
         enddo

! third order
         do l=2,nl
            do r=1,3
            do q=1,3
            do p=1,3
               f(l,i1,j1,k1)=f(l,i1,j1,k1) + h3(p,q,r,l)*a1_3(p,q,r)
            enddo
            enddo
            enddo
         enddo


! scale with weights
!      do l=1,nl
!         f(l,i1,j1,k1)= weights(l)*f(l,i1,j1,k1)
!      enddo
         f( 1,i1,j1,k1)= weights( 1)*f( 1,i1,j1,k1)
         f( 2,i1,j1,k1)= weights( 2)*f( 2,i1,j1,k1)
         f( 3,i1,j1,k1)= weights( 3)*f( 3,i1,j1,k1)
         f( 4,i1,j1,k1)= weights( 4)*f( 4,i1,j1,k1)
         f( 5,i1,j1,k1)= weights( 5)*f( 5,i1,j1,k1)
         f( 6,i1,j1,k1)= weights( 6)*f( 6,i1,j1,k1)
         f( 7,i1,j1,k1)= weights( 7)*f( 7,i1,j1,k1)
         f( 8,i1,j1,k1)= weights( 8)*f( 8,i1,j1,k1)
         f( 9,i1,j1,k1)= weights( 9)*f( 9,i1,j1,k1)
         f(10,i1,j1,k1)= weights(10)*f(10,i1,j1,k1)
         f(11,i1,j1,k1)= weights(11)*f(11,i1,j1,k1)
         f(12,i1,j1,k1)= weights(12)*f(12,i1,j1,k1)
         f(13,i1,j1,k1)= weights(13)*f(13,i1,j1,k1)
         f(14,i1,j1,k1)= weights(14)*f(14,i1,j1,k1)
         f(15,i1,j1,k1)= weights(15)*f(15,i1,j1,k1)
         f(16,i1,j1,k1)= weights(16)*f(16,i1,j1,k1)
         f(17,i1,j1,k1)= weights(17)*f(17,i1,j1,k1)
         f(18,i1,j1,k1)= weights(18)*f(18,i1,j1,k1)
         f(19,i1,j1,k1)= weights(19)*f(19,i1,j1,k1)
#ifndef D3Q19
         f(20,i1,j1,k1)= weights(20)*f(20,i1,j1,k1)
         f(21,i1,j1,k1)= weights(21)*f(21,i1,j1,k1)
         f(22,i1,j1,k1)= weights(22)*f(22,i1,j1,k1)
         f(23,i1,j1,k1)= weights(23)*f(23,i1,j1,k1)
         f(24,i1,j1,k1)= weights(24)*f(24,i1,j1,k1)
         f(25,i1,j1,k1)= weights(25)*f(25,i1,j1,k1)
         f(26,i1,j1,k1)= weights(26)*f(26,i1,j1,k1)
         f(27,i1,j1,k1)= weights(27)*f(27,i1,j1,k1)
#endif
      endif
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
