module m_vreman_kernel
contains
#ifdef _CUDA
   attributes(global)&
#endif
   subroutine vreman_kernel(tau,f, H2, const, kinevisc, eps)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx,ny,nz
   use mod_D3Q27setup, only : nl
   implicit none
   real, intent(in)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(out)   :: tau(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: H2(3,3,nl)
   real, value, intent(in) :: eps
   real, value, intent(in) :: kinevisc
   real, value, intent(in) :: const

   real   :: alpha(3,3)
   real   :: alphamag
   real   :: beta(3,3)
   real   :: Bbeta
   real   :: eddyvisc
   real   :: tmp

   integer :: i, j, k, l, m, q, p

#ifdef _CUDA
   i = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k = threadIdx%z + (blockIdx%z - 1) * blockDim%z
   if (i > nx) return
   if (j > ny) return
   if (k > nz) return
   !eps = sqrt(tiny(1.0))
#else
   !eps = sqrt(tiny(1.0))
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k, l, p, q, m, alpha, alphamag, beta, Bbeta, eddyvisc, tmp)&
!$OMP                           SHARED(f, H2, nl, nx, ny, nz, const, kinevisc, tau, eps)
   do k=1,nz
   do j=1,ny
   do i=1,nx
#endif

! Eq (11) from Jacob 2018 is identical to the 33a from Feng (2021)
      do q=1,3
      do p=1,3
         alpha(p,q) =  0.0
      enddo
      enddo

      do l=1,nl
         tmp=f(l,i,j,k)
         do q=1,3
         do p=1,3
            alpha(p,q) = alpha(p,q) + H2(p,q,l)*tmp
         enddo
         enddo
      enddo

      alphamag=0.00001
      do q=1,3
      do p=1,3
         alphamag=alphamag+alpha(p,q)*alpha(p,q)
      enddo
      enddo

!! beta = del^2 * alpha' * alpha
      do q=1,3
      do p=1,3
         beta(p,q)=0.00001
         do m=1,3
            beta(p,q)=beta(p,q)+alpha(m,p)*alpha(m,q)
         enddo
      enddo
      enddo


      Bbeta=beta(1,1)*beta(2,2) - beta(1,2)**2  &
           +beta(1,1)*beta(3,3) - beta(1,3)**2  &
           +beta(2,2)*beta(3,3) - beta(2,3)**2


      tmp=Bbeta/alphamag
      if (tmp > eps) then
         eddyvisc=const*sqrt(tmp)
      else
         eddyvisc=0.0
      endif
      tau(i,j,k) = 3.0*(kinevisc + eddyvisc) + 0.5

#ifndef _CUDA
    enddo
    enddo
    enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
