module m_vreman
!  Vreman (2004) subgridscale turbulence model
contains

subroutine vreman(f, tau)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile, only : ivreman,kinevisc,smagorinsky
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime
   use m_vreman_kernel
   implicit none
   real, intent(in)      :: f(nl,0:nx+1,0:ny+1,0:nz+1) ! Nonequilibrium f as input
   real, intent(out)     :: tau(0:nx+1,0:ny+1,0:nz+1)              ! Tau including subgrid scale mixing
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: tau
   integer :: tx, ty, tz, bx, by, bz
#endif

   real :: const           ! c in Vreman 2004 Eq (5)
   real :: eps
   integer :: i, j, k

   integer, parameter :: icpu=6
   call cpustart()

   if (ivreman /= 1) then
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

