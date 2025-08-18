module m_vreman
!  Vreman (2004) subgridscale turbulence model
contains

subroutine vreman(f, tau)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile, only : ivreman,kinevisc,smagorinsky
   use m_wtime
   use m_vreman_kernel
   implicit none
   real, intent(in)      :: f(nl,0:nx+1,0:ny+1,0:nz+1) ! Nonequilibrium f as input
   real, intent(out)     :: tau(nx,ny,nz)              ! Tau including subgrid scale mixing
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: tau
#endif

   real :: const           ! c in Vreman 2004 Eq (5)
   real :: tmp
   real :: eps
   integer :: i, j, k, l, m, p, q

   integer, parameter :: icpu=6
   integer :: tx, ty, tz, bx, by, bz
   call cpustart()

   if (ivreman /= 1) then
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k) SHARED(tau, kinevisc)
#endif
      do k=1,nz
         do j=1,ny
            do i=1,nx
               tau(i,j,k) = 3.0*kinevisc + 0.5
            enddo
         enddo
      enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
      return
   endif

   const=2.5*smagorinsky**2
   eps = sqrt(tiny(1.0))

#ifdef _CUDA
   tx=ntx; bx=(nx+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=ntz; bz=(nz+tz-1)/tz
#endif
   call vreman_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(tau, f, H2, const, kinevisc, nx, ny, nz, nl, eps)


   call cpufinish(icpu)

end subroutine
end module

