module m_advection
contains
subroutine advection(tempout,tempin,u,v,w,tau,n)
! f enter routine in feq following collisions
! f is returned in f
   use mod_dimensions, only : nx, ny, nz
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_readinfile, only : iablvisc,istable,p2l
   use m_wtime
   use m_advection_kernel
!   use m_theta_profile_reduce
!   use m_theta_profile_partials
   implicit none
   integer, intent(in) :: n
   real, intent(inout) :: tempin (n,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: tempout(n,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: u(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: v(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: w(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: tau(0:nx+1,0:ny+1,0:nz+1)

   real            :: weights(-1:1,-1:1,-1:1)
   real, parameter :: weights_h(-1:1,-1:1,-1:1) = reshape(&
                         [ 0.125, 0.250, 0.125, 0.250, 1.000, 0.250, 0.125, 0.250, 0.125, &
                           0.250, 1.000, 0.250, 1.000,-10.000, 1.000, 0.250, 1.000, 0.250, &
                           0.125, 0.250, 0.125, 0.250, 1.000, 0.250, 0.125, 0.250, 0.125 ], shape(weights_h) )
#ifdef _CUDA
   attributes(device) :: tempin
   attributes(device) :: tempout
   attributes(device) :: u,v,w,tau
   attributes(device) :: weights
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer, parameter :: icpu=30
   integer i,j,k
   integer :: isel1, isel2, jsel1, jsel2
   integer :: nsel
   call cpustart()

   weights=weights_h

#ifdef _CUDA
   tx=ntx; bx=(n*(nx+1)+tx-1)/tx
   ty=nty; by=(ny+ty-1)/ty
   tz=ntz; bz=(nz+tz-1)/tz
#endif
   call advection_kernel&
#ifdef _CUDA
        &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
        &(tempin, tempout, u, v, w, weights, tau, n)

   call cpufinish(icpu)

end subroutine
end module
