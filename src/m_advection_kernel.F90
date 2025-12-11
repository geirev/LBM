module m_advection_kernel
contains

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine advection_kernel(tracerin, tracerout, u, v, w, weights, tau, n)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx, ny, nz
   implicit none

   ! f and feq flattened in (l,i), explicit dims in j,k
   integer, value    :: n
   real, intent(in)  :: tracerin (n*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(out) :: tracerout(n*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(in)  :: u(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: v(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: w(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: tau(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: weights(-1:1,-1:1,-1:1)

   real dtx,dty,dtz,diff,uu,vv,ww,tt
   integer :: idx, j, k, i, l, ii, jj, kk
   integer :: idxp, idxm
   real, parameter :: kappa_mol = 1.0e-3

#ifdef _CUDA
   idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x   ! flattened li index
   j   = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k   = threadIdx%z + (blockIdx%z - 1) * blockDim%z

   if (idx < n+1 .or. idx > n*(nx+1)) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(idx,i,l,j,k,dtx,dty,dtz,diff,uu,vv,ww,tt,idxm,idxp) &
!$OMP                        &SHARED(u,v,w,tracerin,tracerout,tau,weights)
   do k = 1, nz
   do j = 1, ny
   do idx = n+1, n*(nx+1)
#endif
      l  = mod(idx-1, n) + 1
      i  = (idx-1)/n

      uu=u(i,j,k)
      vv=v(i,j,k)
      ww=w(i,j,k)

      tt = (tau(i,j,k) - 0.5)/3.0
      tt = max(tt, 0.0)
      tt = min(tt, 0.16)
      tt = tt + kappa_mol


!--------------------------------------------------------------
! (1) Predictor stage: θ* = θ^n + dt * RHS(θ^n)
!--------------------------------------------------------------
! dx = dy = dz = dt = 1 in lattice units, so dt/(dx) = 1 etc.
! RHS advection with centered differences
!!      dtx=tracerin(idx+n,j  ,k  )   - tracerin(idx-n ,j  ,k  )
!!      dty=tracerin(idx  ,j+1,k  )   - tracerin(idx   ,j-1,k  )
!!      dtz=tracerin(idx  ,j  ,k+1)   - tracerin(idx   ,j  ,k-1)

! Advection with simple 1st-order upwind

! x-direction: idx +/- n moves one i-cell (since l is fastest)
   idxp = idx + n
   idxm = idx - n

   if (uu > 0.0) then
      dtx = tracerin(idx  ,j,k) - tracerin(idxm,j,k)
   else
      dtx = tracerin(idxp,j,k) - tracerin(idx  ,j,k)
   end if

   ! y-direction
   if (vv > 0.0) then
      dty = tracerin(idx,j  ,k) - tracerin(idx,j-1,k)
   else
      dty = tracerin(idx,j+1,k) - tracerin(idx,j  ,k)
   end if

   ! z-direction
   if (ww > 0.0) then
      dtz = tracerin(idx,j,k  ) - tracerin(idx,j,k-1)
   else
      dtz = tracerin(idx,j,k+1) - tracerin(idx,j,k  )
   end if

! RHS diffusion
      diff=0.0
      do kk=-1,1
      do jj=-1,1
      do ii=-1,1
         diff = diff + weights(ii,jj,kk)*tracerin(idx+ii*n,j+jj,k+kk)
      end do
      end do
      end do

! predictor centered
!      tracerout(idx,j,k) = tracerin(idx,j,k) -0.5*(uu*dtx + vv*dty + ww*dtz) + tt*diff
! predictor upwind
      tracerout(idx,j,k) = tracerin(idx,j,k) -(uu*dtx + vv*dty + ww*dtz) + tt*diff


!--------------------------------------------------------------
! (2) Corrector stage: θ^{n+1} = 0.5θ^n + 0.5(θ* + dt*RHS(θ*))
!--------------------------------------------------------------

! RHS advection
!!      dtx=tracerout(idx+n,j  ,k  )   - tracerout(idx-n,j  ,k  )
!!      dty=tracerout(idx  ,j+1,k  )   - tracerout(idx  ,j-1,k  )
!!      dtz=tracerout(idx  ,j  ,k+1)   - tracerout(idx  ,j  ,k-1)

! RHS advection (corrector stage) - upwind again
      if (uu > 0.0) then
         dtx = tracerout(idx  ,j,k) - tracerout(idxm,j,k)
      else
         dtx = tracerout(idxp,j,k) - tracerout(idx  ,j,k)
      end if

      if (vv > 0.0) then
         dty = tracerout(idx,j  ,k) - tracerout(idx,j-1,k)
      else
         dty = tracerout(idx,j+1,k) - tracerout(idx,j  ,k)
      end if

      if (ww > 0.0) then
         dtz = tracerout(idx,j,k  ) - tracerout(idx,j,k-1)
      else
         dtz = tracerout(idx,j,k+1) - tracerout(idx,j,k  )
      end if

! RHS diffusion
      diff=0.0
      do kk=-1,1
      do jj=-1,1
      do ii=-1,1
         diff = diff + weights(ii,jj,kk)*tracerout(idx+ii*n,j+jj,k+kk)
      end do
      end do
      end do

! Corrector step centered
!      tracerout(idx,j,k) = 0.5*tracerin(idx,j,k) + 0.5*( tracerout(idx,j,k) - 0.5*(uu*dtx+vv*dty+ww*dtz) + tt*diff)
! Corrector step upwind
      tracerout(idx,j,k) = 0.5*tracerin(idx,j,k) + 0.5*( tracerout(idx,j,k) - (uu*dtx+vv*dty+ww*dtz) + tt*diff )

#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
