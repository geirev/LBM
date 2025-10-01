module m_postcoll
contains

subroutine postcoll(f, feq, tau, rho, u, v, w)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_readinfile, only : ibgk,ihrr
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz, ivreman,kinevisc,smagorinsky
#endif
   use m_wtime
   use m_postcoll_kernel

   implicit none
   real, intent(in)      :: rho(nx,ny,nz)
   real, intent(in)      :: u(nx,ny,nz)
   real, intent(in)      :: v(nx,ny,nz)
   real, intent(in)      :: w(nx,ny,nz)
   real, intent(out)     :: feq(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)   :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout)   :: tau(0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: feq
   attributes(device) :: f
   attributes(device) :: tau
   integer :: tx, ty, tz, bx, by, bz
#endif

   real :: const           ! c in Vreman 2004 Eq (5)
   real :: eps

   real, parameter :: inv1cs2 = 1.0/(cs2)
   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv2cs6 = 1.0/(2.0*cs6)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
   integer, parameter :: icpu=18

   const=2.5*smagorinsky**2
   eps = sqrt(tiny(1.0))

   call cpustart()
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call postcoll_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(f, feq, tau, rho, u, v, w, inv1cs2, inv2cs4, inv6cs6, eps, kinevisc, const, ibgk, ihrr, ivreman)

   call cpufinish(icpu)

end subroutine
end module
