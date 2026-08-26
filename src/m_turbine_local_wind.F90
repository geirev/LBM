module m_turbine_local_wind
contains

subroutine turbine_local_wind(turbine,u,v,w,rho,dist_upstream_D,nsrad, &
                              uavg,vavg,wavg,speed,winddir,unormal)
! Compute the disk-averaged local wind upstream of one turbine.
!
! The sampling plane:
!  - is parallel to the current rotor plane,
!  - has diameter equal to the rotor diameter,
!  - is centered dist_upstream_D rotor diameters upstream,
!  - uses the current turbine yaw and tilt.
!
! The returned direction is
!
!       winddir = atan2(vavg,uavg)
!
! and is therefore in radians, consistent with turbine yaw.
!
! Recommended initial values:
!
!       dist_upstream_D = 1.0
!       nsrad            = 6
!
! nsrad=6 gives approximately pi*6^2 = 113 sampling points.

   use mod_dimensions
   use mod_turbine_def, only : radiusnd
   use mod_turbines, only : turbine_t
   use m_turbine_rotor_basis
   use m_turbine_local_wind_kernel
#ifdef _CUDA
   use cudafor
#endif
#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : j_start,j_end
#endif

   implicit none

#ifndef MPI
   integer, parameter :: j_start = 1
   integer, parameter :: j_end   = ny
#endif

   type(turbine_t), intent(in) :: turbine

   real, intent(in) :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: w  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)

#ifdef _CUDA
   attributes(device) :: u,v,w,rho
#endif

   real,    intent(in) :: dist_upstream_D
   integer, intent(in) :: nsrad

   real, intent(out) :: uavg,vavg,wavg
   real, intent(out) :: speed
   real, intent(out) :: winddir
   real, intent(out) :: unormal

   integer :: nside,npt
#ifdef _CUDA
   integer :: tpb,nblocks,istat
#endif
#ifdef MPI
   integer :: ierr
#endif

   real :: diameter
   real :: e_axis(3),e1(3),e2(3)
   real :: xc,yc,zc

   real :: sums(4)

#ifdef _CUDA
   real, device :: sums_d(4)
#endif


!------------------------------------------------------------------
! Check sampling resolution.
!------------------------------------------------------------------
   if (nsrad < 1) then
      uavg   = 0.0
      vavg   = 0.0
      wavg   = 0.0
      speed  = 0.0
      winddir = turbine%yaw
      return
   endif


!------------------------------------------------------------------
! Rotor basis corresponding to the current turbine yaw and tilt.
!------------------------------------------------------------------
   call turbine_rotor_basis(turbine%yaw,turbine%tilt,e_axis,e1,e2)

   diameter = 2.0*radiusnd


!------------------------------------------------------------------
! Center of the upstream sampling disk.
!
! Positive e_axis points along the turbine rotor axis, so the
! measurement disk is placed upstream along -e_axis.
!------------------------------------------------------------------
   xc = turbine%xhub - dist_upstream_D*diameter*e_axis(1)
   yc = turbine%yhub - dist_upstream_D*diameter*e_axis(2)
   zc = turbine%zhub - dist_upstream_D*diameter*e_axis(3)


!------------------------------------------------------------------
! Number of candidate points on the square grid surrounding the
! circular sampling disk.
!------------------------------------------------------------------
   nside = 2*nsrad + 1
   npt   = nside*nside


!------------------------------------------------------------------
! Accumulate local velocity sums.
!------------------------------------------------------------------
   sums = 0.0

#ifdef _CUDA
   sums_d = 0.0

   tpb     = 128
   nblocks = (npt + tpb - 1)/tpb

   call turbine_local_wind_kernel&
#ifdef _CUDA
   &<<<nblocks,tpb>>>&
#endif
   &(u,v,w,rho, &
     xc,yc,zc, &
     e1(1),e1(2),e1(3), &
     e2(1),e2(2),e2(3), &
     radiusnd,nsrad,j_start,j_end,sums_d)

   istat = cudaDeviceSynchronize()

   sums = sums_d

#else

   call turbine_local_wind_kernel&
   &(u,v,w,rho, &
     xc,yc,zc, &
     e1(1),e1(2),e1(3), &
     e2(1),e2(2),e2(3), &
     radiusnd,nsrad,j_start,j_end,sums)

#endif


!------------------------------------------------------------------
! Sum contributions across MPI tiles.
!
! Each sampling point is evaluated on exactly one rank, determined
! from the global lower j-index of its interpolation stencil.
!------------------------------------------------------------------
#ifdef MPI
   call MPI_Allreduce(MPI_IN_PLACE,sums,4,MPI_REAL,MPI_SUM, &
                      MPI_COMM_WORLD,ierr)
#endif


!------------------------------------------------------------------
! Disk-averaged velocity vector.
!------------------------------------------------------------------
   if (sums(4) > 0.0) then

      uavg = sums(1)/sums(4)
      vavg = sums(2)/sums(4)
      wavg = sums(3)/sums(4)

   else

      uavg = 0.0
      vavg = 0.0
      wavg = 0.0

   endif


!------------------------------------------------------------------
! Magnitude of the averaged velocity vector.
!------------------------------------------------------------------
   speed = sqrt(uavg*uavg + vavg*vavg + wavg*wavg)


!------------------------------------------------------------------
! Horizontal wind direction in radians.
!
! Averaging u and v before atan2 avoids all 0/360-degree wrapping
! problems associated with averaging direction angles directly.
!------------------------------------------------------------------
   if (uavg /= 0.0 .or. vavg /= 0.0) then
      winddir = atan2(vavg,uavg)
   else
      winddir = turbine%yaw
   endif

   unormal = uavg*e_axis(1) + &
             vavg*e_axis(2) + &
             wavg*e_axis(3)

end subroutine
end module
