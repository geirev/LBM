module m_vreman_tau_kernel
! Eq (11) from Jacob 2018 is identical to the 33a from Feng (2021)
#ifdef _CUDA
   use cudafor
#endif
   implicit none
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine vreman_tau_kernel(tau, eddyvisc, Bbeta, alphamag, kinevisc, const, nx, ny, nz)
   implicit none
   integer, value      :: nx, ny, nz
   real, intent(out)   :: tau(nx,ny,nz)
   real, intent(out)   :: eddyvisc(nx,ny,nz)
   real, intent(in)    :: Bbeta(nx,ny,nz)
   real, intent(in)    :: alphamag(nx,ny,nz)
   real, value         :: kinevisc
   real, value         :: const
   integer :: i, j, k
#ifdef _CUDA
   attributes(device) :: tau
   attributes(device) :: eddyvisc
   attributes(device) :: Bbeta
   attributes(device) :: alphamag
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx) return
   if (j > ny) return
   if (k > nz) return
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, p, q, m) SHARED(tau, eddyvisc, Bbeta, alphamag, kinevisc, const)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif
      eddyvisc(i,j,k)=const*sqrt(Bbeta(i,j,k)/alphamag(i,j,k))
      tau(i,j,k) = 3.0*(kinevisc + eddyvisc(i,j,k)) + 0.5
#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module

