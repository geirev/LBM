! Equilibrium distribution \citet{fen21a} Eq. (32) or jac18a eq (27)
module m_fequil3_1ord_kernel
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine fequil3_1ord_kernel(feq, rho, vel, nx, ny, nz, nl, cxs, cys, czs, cs2, inv1cs2)
   implicit none
   integer, value       :: nx, ny, nz, nl
   real, intent(out)    :: feq(nl,nx+2,ny+2,nz+2)
   real, intent(in)     :: rho(nx,ny,nz)
   real, intent(in)     :: vel(3,nx,ny,nz)
   integer, intent(in)  :: cxs(nl)
   integer, intent(in)  :: cys(nl)
   integer, intent(in)  :: czs(nl)
   real, value          :: inv1cs2
   real, value          :: cs2
   integer :: i, j, k, l
#ifdef _CUDA
   attributes(device) :: feq
   attributes(device) :: rho
   attributes(device) :: vel
   attributes(device) :: cxs
   attributes(device) :: cys
   attributes(device) :: czs
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx .or. j > ny .or. k > nz) return
#else
!$OMP PARALLEL DO collapse(3) DEFAULT(NONE) PRIVATE(i, j, k) SHARED(f, vel, rho, nx, ny, nz, cs2, inv1cs2)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      do l=1,nl
         feq(l,i+1,j+1,k+1) = rho(i,j,k) * (cs2 + real(cxs(l))*vel(1,i,j,k)&
                                                + real(cys(l))*vel(2,i,j,k)&
                                                + real(czs(l))*vel(3,i,j,k)) * inv1cs2
      enddo

!      feq( 1,i+1,j+1,k+1) = rho(i,j,k) * (cs2                                                ) * inv1cs2
!      feq( 2,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)                                ) * inv1cs2
!      feq( 3,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)                                ) * inv1cs2
!      feq( 4,i+1,j+1,k+1) = rho(i,j,k) * (cs2                  + vel(2,i,j,k)                ) * inv1cs2
!      feq( 5,i+1,j+1,k+1) = rho(i,j,k) * (cs2                  - vel(2,i,j,k)                ) * inv1cs2
!      feq( 6,i+1,j+1,k+1) = rho(i,j,k) * (cs2                                  - vel(3,i,j,k)) * inv1cs2
!      feq( 7,i+1,j+1,k+1) = rho(i,j,k) * (cs2                                  + vel(3,i,j,k)) * inv1cs2
!      feq( 8,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  + vel(2,i,j,k)                ) * inv1cs2
!      feq( 9,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  - vel(2,i,j,k)                ) * inv1cs2
!      feq(10,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  - vel(2,i,j,k)                ) * inv1cs2
!      feq(11,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  + vel(2,i,j,k)                ) * inv1cs2
!      feq(12,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)                  - vel(3,i,j,k)) * inv1cs2
!      feq(13,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)                  + vel(3,i,j,k)) * inv1cs2
!      feq(14,i+1,j+1,k+1) = rho(i,j,k) * (cs2                  + vel(2,i,j,k)  + vel(3,i,j,k)) * inv1cs2
!      feq(15,i+1,j+1,k+1) = rho(i,j,k) * (cs2                  - vel(2,i,j,k)  - vel(3,i,j,k)) * inv1cs2
!      feq(16,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)                  + vel(3,i,j,k)) * inv1cs2
!      feq(17,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)                  - vel(3,i,j,k)) * inv1cs2
!      feq(18,i+1,j+1,k+1) = rho(i,j,k) * (cs2                  - vel(2,i,j,k)  + vel(3,i,j,k)) * inv1cs2
!      feq(19,i+1,j+1,k+1) = rho(i,j,k) * (cs2                  + vel(2,i,j,k)  - vel(3,i,j,k)) * inv1cs2
!      feq(20,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  + vel(2,i,j,k)  + vel(3,i,j,k)) * inv1cs2
!      feq(21,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  - vel(2,i,j,k)  - vel(3,i,j,k)) * inv1cs2
!      feq(22,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  - vel(2,i,j,k)  - vel(3,i,j,k)) * inv1cs2
!      feq(23,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  + vel(2,i,j,k)  + vel(3,i,j,k)) * inv1cs2
!      feq(24,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  + vel(2,i,j,k)  - vel(3,i,j,k)) * inv1cs2
!      feq(25,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  - vel(2,i,j,k)  + vel(3,i,j,k)) * inv1cs2
!      feq(26,i+1,j+1,k+1) = rho(i,j,k) * (cs2  - vel(1,i,j,k)  + vel(2,i,j,k)  - vel(3,i,j,k)) * inv1cs2
!      feq(27,i+1,j+1,k+1) = rho(i,j,k) * (cs2  + vel(1,i,j,k)  - vel(2,i,j,k)  + vel(3,i,j,k)) * inv1cs2
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
