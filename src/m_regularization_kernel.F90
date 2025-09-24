module m_regularization_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine regularization_kernel(f, feq, u, v, w, nx, ny, nz, nl, h2, h3, weights, inv2cs4, inv6cs6)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   integer, value      :: nx, ny, nz, nl
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)
   real, intent(in)    :: h2(3,3,nl)
   real, intent(in)    :: h3(3,3,3,nl)
   real, intent(in)    :: weights(nl)
   real, value         :: inv2cs4
   real, value         :: inv6cs6


   real :: vel(3)
   real :: a1_2(3,3)
   real :: a1_3(3,3,3)
   real :: tmp

   integer :: i, j, k, l, p, q, r, i1, j1, k1
#ifdef _CUDA
   i = threadidx%x + (blockidx%x - 1) * blockdim%x
   j = threadidx%y + (blockidx%y - 1) * blockdim%y
   k = threadidx%z + (blockidx%z - 1) * blockdim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(3) DEFAULT(none) PRIVATE(i, j, k, l, p, q, r, i1, j1, k1, vel, a1_2, a1_3, tmp)&
!$OMP             & SHARED(f, u, v, w, nx, ny, nz, nl, h2, h3, weights, inv2cs4, inv6cs6)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
! copy u,v,w to vel(1:3)
      vel(1)=u(i,j,k)
      vel(2)=v(i,j,k)
      vel(3)=w(i,j,k)

! computing a1_2
      a1_2(:,:)=0.0
      do l=1,nl
         tmp=f(l,i,j,k)-feq(l,i,j,k)
         do q=1,3
         do p=1,3
            a1_2(p,q) = a1_2(p,q) + h2(p,q,l)*tmp
         enddo
         enddo
      enddo


! computing a1_3 using Maspalinas (2015, appedix B)
      do r=1,3
      do q=1,3
      do p=1,3
         a1_3(p,q,r)=vel(p)*a1_2(q,r) + vel(q)*a1_2(r,p) +  vel(r)*a1_2(p,q)
      enddo
      enddo
      enddo

! compute actual 3rd-order Hermite coefficients: a1_3(p,q,r) = sum_l h3(p,q,r,l) * f(l)
!        a1_3(:,:,:) = 0.0
!        do l = 1, nl
!           do r = 1, 3
!           do q = 1, 3
!           do p = 1, 3
!              a1_3(p,q,r) = a1_3(p,q,r) + h3(p,q,r,l) * f(l,i1,j1,k1)
!           end do
!           end do
!           end do
!        end do


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

      do l=1,nl
         tmp = 0.0
         do q=1,3
         do p=1,3
            tmp = tmp + h2(p,q,l)*a1_2(p,q)
         end do
         end do

         do r=1,3
         do q=1,3
         do p=1,3
            tmp = tmp + h3(p,q,r,l)*a1_3(p,q,r)
         end do
         end do
         end do

         f(l,i,j,k) = weights(l) * tmp
      end do
#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
