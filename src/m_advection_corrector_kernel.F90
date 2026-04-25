module m_advection_corrector_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
   subroutine advection_corrector_kernel(tracerin, tracerstar, tracerout, u, v, w, weights, tau, n)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx, ny, nz
   implicit none

   integer, value    :: n
   real, intent(in)  :: tracerin  (n*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(in)  :: tracerstar(n*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(out) :: tracerout (n*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(in)  :: u(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: v(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: w(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: tau(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: weights(-1:1,-1:1,-1:1)

   integer :: idx, i, l, j, k, ii, jj, kk, idxp, idxm
   real :: uu, vv, ww, tt, dtx, dty, dtz, diff
   real, parameter :: kappa_mol = 1.0e-3

#ifdef _CUDA
   idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x
   j   = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k   = threadIdx%z + (blockIdx%z - 1) * blockDim%z

   if (idx < n+1 .or. idx > n*(nx+1)) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(idx,i,l,j,k,ii,jj,kk,idxp,idxm,uu,vv,ww,tt,dtx,dty,dtz,diff) &
!$OMP SHARED(tracerin,tracerstar,tracerout,u,v,w,tau,weights)
   do k = 1, nz
   do j = 1, ny
   do idx = n+1, n*(nx+1)
#endif

      l = mod(idx-1, n) + 1
      i = (idx-1)/n

      uu = u(i,j,k)
      vv = v(i,j,k)
      ww = w(i,j,k)

      tt = (tau(i,j,k) - 0.5)/3.0
      tt = max(tt, 0.0)
      tt = min(tt, 0.16)
      tt = tt + kappa_mol

      idxp = idx + n
      idxm = idx - n

      if (uu > 0.0) then
         dtx = tracerstar(idx,j,k) - tracerstar(idxm,j,k)
      else
         dtx = tracerstar(idxp,j,k) - tracerstar(idx,j,k)
      end if

      if (vv > 0.0) then
         dty = tracerstar(idx,j,k) - tracerstar(idx,j-1,k)
      else
         dty = tracerstar(idx,j+1,k) - tracerstar(idx,j,k)
      end if

      if (ww > 0.0) then
         dtz = tracerstar(idx,j,k) - tracerstar(idx,j,k-1)
      else
         dtz = tracerstar(idx,j,k+1) - tracerstar(idx,j,k)
      end if

      diff = 0.0
      do kk = -1, 1
      do jj = -1, 1
      do ii = -1, 1
         diff = diff + weights(ii,jj,kk) * tracerstar(idx+ii*n,j+jj,k+kk)
      end do
      end do
      end do

      tracerout(idx,j,k) = 0.5*tracerin(idx,j,k)                       &
                         + 0.5*( tracerstar(idx,j,k)                   &
                         - (uu*dtx + vv*dty + ww*dtz)                  &
                         + tt*diff )

#ifndef _CUDA
   end do
   end do
   end do
!$OMP END PARALLEL DO
#endif

   end subroutine advection_corrector_kernel

end module m_advection_corrector_kernel
