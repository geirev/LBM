module m_turbine_local_wind_kernel
contains
#ifdef _CUDA
   attributes(global) &
#endif
subroutine turbine_local_wind_kernel(u,v,w,rho, &
                                     xc,yc,zc, &
                                     e1x,e1y,e1z, &
                                     e2x,e2y,e2z, &
                                     radius,nsrad,j_start,j_end,sums)
! Compute contributions to the disk-averaged upstream wind velocity.
!
! The sampling disk:
!  - is centered at (xc,yc,zc),
!  - is parallel to the rotor plane,
!  - has radius "radius",
!  - is discretized by a Cartesian grid in the local e1-e2 plane.
!
! sums(1) = sum(u)
! sums(2) = sum(v)
! sums(3) = sum(w)
! sums(4) = number of valid sampling points
!
! Under MPI, each sampling point is evaluated only by the rank owning
! the lower global j-index required for interpolation.

#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions
   use m_turbine_interpolate_velocity

   implicit none

   real, intent(in) :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: w  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)

   real, value    :: xc,yc,zc
   real, value    :: e1x,e1y,e1z
   real, value    :: e2x,e2y,e2z
   real, value    :: radius
   integer, value :: nsrad
   integer, value :: j_start,j_end

   real, intent(inout) :: sums(4)

   integer :: p,nside,npt
   integer :: i1,i2
   integer :: ig,jg,kg

   real :: ds,s1,s2
   real :: xg,yg,zg
   real :: ux,uy,uz,dens
#ifdef _CUDA
   real :: old
#endif


   nside = 2*nsrad + 1
   npt   = nside*nside

#ifdef _CUDA
   p = (blockIdx%x-1)*blockDim%x + threadIdx%x
   if (p < 1 .or. p > npt) return
#else
   do p = 1,npt
#endif

      !---------------------------------------------------------------
      ! Convert the one-dimensional point number into coordinates
      ! on the square sampling grid.
      !---------------------------------------------------------------
      i1 = mod(p-1,nside) - nsrad
      i2 = (p-1)/nside    - nsrad

      ds = radius/real(nsrad)

      s1 = real(i1)*ds
      s2 = real(i2)*ds

      !---------------------------------------------------------------
      ! Retain only points inside the circular rotor disk.
      !---------------------------------------------------------------
      if (s1*s1 + s2*s2 <= radius*radius) then

         !------------------------------------------------------------
         ! Sampling-point position in global lattice coordinates:
         !
         !    x = center + s1*e1 + s2*e2
         !------------------------------------------------------------
         xg = xc + s1*e1x + s2*e2x
         yg = yc + s1*e1y + s2*e2y
         zg = zc + s1*e1z + s2*e2z

         ig = floor(xg)
         jg = floor(yg)
         kg = floor(zg)

         !------------------------------------------------------------
         ! Sampling point must have a valid interpolation stencil.
         !
         ! In i and k, the upper interpolation point may use the
         ! nx+1/nz+1 ghost cell.
         !
         ! In j, ownership is determined using the global lower cell.
         !------------------------------------------------------------
         if (ig >= 1 .and. ig <= nx .and. &
             kg >= 1 .and. kg <= nz .and. &
             jg >= j_start .and. jg <= j_end) then

            call turbine_interpolate_velocity(u,v,w,rho, &
                 xg,yg,zg,j_start,ux,uy,uz,dens)

#ifdef _CUDA
            old = atomicadd(sums(1),ux)
            old = atomicadd(sums(2),uy)
            old = atomicadd(sums(3),uz)
            old = atomicadd(sums(4),1.0)
#else
            sums(1) = sums(1) + ux
            sums(2) = sums(2) + uy
            sums(3) = sums(3) + uz
            sums(4) = sums(4) + 1.0
#endif

         endif
      endif

#ifndef _CUDA
   enddo
#endif

end subroutine
end module
