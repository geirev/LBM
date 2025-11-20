module m_turbine_point_forces_gpu
!    GPU wrapper: allocate device arrays for points and Fvec,
!    launch CUDA kernel over points, then copy forces back.
contains
subroutine turbine_point_forces_gpu(points_global, rho, u, v, w, Fvec_local, np)
   use mod_dimensions, only : nx, ny, nz
   use mod_turbines,   only : point_t
   use m_turbine_point_forces_kernel
#ifdef _CUDA
   use cudafor
#endif
#ifdef MPI
   use m_mpi_decomp_init, only : j_start, j_end
#endif
   implicit none
#ifndef MPI
   integer, parameter :: j_start = 1, j_end = ny
#endif

   type(point_t), intent(in) :: points_global(:)
   real, intent(in)          :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)          :: w  (0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device)        :: u,v,w,rho
#endif
   integer, intent(in)       :: np

   real, intent(out)         :: Fvec_local(3, np)

   type(point_t), allocatable :: points(:)
   real,          allocatable :: Fvec(:,:)
#ifdef _CUDA
   attributes(device) :: points
   attributes(device) :: Fvec
#endif

   integer :: tpb, nblocks, istat

   if (np <= 0) then
      Fvec_local = 0.0
      return
   end if

   allocate(points(np))
   points = points_global

   allocate(Fvec(3,np))
   Fvec = 0.0

   tpb     = 128
   nblocks = (np + tpb - 1) / tpb

   call turbine_point_forces_kernel&
#ifdef _CUDA
        &<<<nblocks, tpb>>>&
#endif
        &(points, np, rho, u, v, w, j_start, j_end, Fvec)

#ifdef _CUDA
   istat = cudaDeviceSynchronize()
#endif

   Fvec_local = Fvec

   deallocate(points, Fvec)
end subroutine
end module
