module m_advection_kernel
contains

#ifdef _CUDA
   attributes(global) &
#endif
   subroutine advection_kernel(tracerin, tracerout, u, v, w, weights, tau)
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx, ny, nz, ntracer
   implicit none

   ! f and feq flattened in (l,i), explicit dims in j,k
   real, intent(in)  :: tracerin (ntracer*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(out) :: tracerout(ntracer*(nx+2), 0:ny+1, 0:nz+1)
   real, intent(in)  :: u(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: v(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: w(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: tau(0:nx+1, 0:ny+1, 0:nz+1)
   real, intent(in)  :: weights(-1:1,-1:1,-1:1)

   real dtx,dty,dtz,diff,uu,vv,ww,tt
   integer :: idx, j, k, i, l, ii, jj, kk

#ifdef _CUDA
   idx = threadIdx%x + (blockIdx%x - 1) * blockDim%x   ! flattened li index
   j   = threadIdx%y + (blockIdx%y - 1) * blockDim%y
   k   = threadIdx%z + (blockIdx%z - 1) * blockDim%z

   if (idx .le. ntracer .or. idx > ntracer*(nx+1)) return
   if (j < 1 .or. j > ny) return
   if (k < 1 .or. k > nz) return
#else
!$OMP PARALLEL DO COLLAPSE(3) PRIVATE(idx,i,l,j,k,dtx,dty,dtz,diff,uu,vv,ww,tt) SHARED(u,v,w,tracerin.tracerout,tau,weights)
   do k = 1, nz
   do j = 1, ny
   do idx = ntracer+1, ntracer*(nx+1)
#endif
      l  = mod(idx-1, ntracer) + 1
      i  = (idx-1)/ntracer

      uu=u(i,j,k)
      vv=v(i,j,k)
      ww=w(i,j,k)
      tt=(tau(i,j,k)-0.5)/3.0

!--------------------------------------------------------------
! (1) Predictor stage: θ* = θ^n + dt * RHS(θ^n)
!--------------------------------------------------------------
! RHS advection
      dtx=tracerin(idx+ntracer,j  ,k  )   - tracerin(idx-ntracer,j  ,k  )
      dty=tracerin(idx        ,j+1,k  )   - tracerin(idx        ,j-1,k  )
      dtz=tracerin(idx        ,j  ,k+1)   - tracerin(idx        ,j  ,k-1)

! RHS diffusion
      diff=0.0
      do kk=-1,1
      do jj=-1,1
      do ii=-1,1
         diff = diff + weights(ii,jj,kk)*tracerin(idx+ii*ntracer,j+jj,k+kk)
      end do
      end do
      end do

! predictor
      tracerout(idx,j,k) = tracerin(idx,j,k) -0.5*(uu*dtx + vv*dty + ww*dtz) + tt*diff


!--------------------------------------------------------------
! (2) Corrector stage: θ^{n+1} = 0.5θ^n + 0.5(θ* + dt*RHS(θ*))
!--------------------------------------------------------------

! RHS advection
      dtx=tracerout(idx+ntracer,j  ,k  )   - tracerout(idx-ntracer,j  ,k  )
      dty=tracerout(idx        ,j+1,k  )   - tracerout(idx        ,j-1,k  )
      dtz=tracerout(idx        ,j  ,k+1)   - tracerout(idx        ,j  ,k-1)

! RHS diffusion
      diff=0.0
      do kk=-1,1
      do jj=-1,1
      do ii=-1,1
         diff = diff + weights(ii,jj,kk)*tracerout(idx+ii*ntracer,j+jj,k+kk)
      end do
      end do
      end do

! Corrector step
      tracerout(idx,j,k) = 0.5*tracerin(idx,j,k) + 0.5*( tracerout(idx,j,k) - 0.5*(uu*dtx+vv*dty+ww*dtz) + tt*diff)

#ifndef _CUDA
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
#endif

end subroutine
end module
