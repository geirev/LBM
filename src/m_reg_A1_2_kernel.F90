module m_reg_A1_2_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine reg_A1_2_kernel(A1_2, H2, f, nx, ny, nz, nl)
   implicit none
   integer, value :: nx, ny, nz, nl
   real, intent(out) :: A1_2(3,3,nx,ny,nz)
   real, intent(in)  :: H2(3,3,nl)
   real, intent(in)  :: f(nl, nx+2, ny+2, nz+2)
   integer :: i, j, k, p, q, l
#ifdef _CUDA
   attributes(device) :: A1_2
   attributes(device) :: H2
   attributes(device) :: f
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO collapse(1) DEFAULT(NONE) PRIVATE(i, j, k, l, p, q) SHARED(f, H2, A1_2, nx, ny, nz, nl)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      l=1
      do q=1,3
      do p=1,3
         A1_2(p,q,i,j,k) = H2(p,q,l)*f(l,i+1,j+1,k+1)
      enddo
      enddo

      do l=2,nl
         do q=1,3
         do p=1,3
            A1_2(p,q,i,j,k) = A1_2(p,q,i,j,k) + H2(p,q,l)*f(l,i+1,j+1,k+1)
         enddo
         enddo
      enddo
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif


!! #ifdef _CUDA
!! !$cuf kernel do(2) <<<*,*>>>
!! #else
!! !$OMP PARALLEL DO collapse(1) DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, r) SHARED(f, H2, A1_2, A1_3)
!! #endif
!!       do k=1,nz
!!       do j=1,ny
!!       do i=1,nx
!!          l=1
!!          do q=1,3
!!          do p=1,3
!!             A1_2(p,q,i,j,k) = H2(p,q,l)*f(l,i,j,k)
!!          enddo
!!          enddo
!!
!!          do l=2,nl
!!             do q=1,3
!!             do p=1,3
!!                A1_2(p,q,i,j,k) = A1_2(p,q,i,j,k) + H2(p,q,l)*f(l,i,j,k)
!!             enddo
!!             enddo
!!          enddo
!!       enddo
!!       enddo
!!       enddo
!! #ifndef _CUDA
!! !$OMP END PARALLEL DO
!! #endif
end subroutine
end module
