module m_turbine_point_kinematics_gpu
contains
subroutine turbine_point_kinematics_gpu(points_global,rho,u,v,w,kin_local,np)
   use mod_dimensions, only : nx,ny,nz
   use mod_turbines, only : point_t
   use m_turbine_point_kinematics_kernel, only : turbine_point_kinematics_kernel
#ifdef _CUDA
   use cudafor
#endif
#ifdef MPI
   use m_mpi_decomp_init, only : j_start,j_end
#endif
   implicit none
#ifndef MPI
   integer, parameter :: j_start=1,j_end=ny
#endif

   type(point_t), intent(in) :: points_global(:)
   real, intent(in) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: u  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: v  (0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: w  (0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: rho,u,v,w
#endif
   integer, intent(in) :: np
   real, intent(out) :: kin_local(3,np)

   type(point_t), allocatable :: points(:)
   real, allocatable :: kin(:,:)
#ifdef _CUDA
   attributes(device) :: points,kin
#endif
   integer :: tpb,nblocks,istat

   if (np<=0) then
      kin_local=0.0
      return
   endif

   allocate(points(np),kin(3,np))
   points=points_global
   kin=0.0
   tpb=128
   nblocks=(np+tpb-1)/tpb

   call turbine_point_kinematics_kernel&
#ifdef _CUDA
        &<<<nblocks,tpb>>>&
#endif
        &(points,np,rho,u,v,w,j_start,j_end,kin)
#ifdef _CUDA
   istat=cudaDeviceSynchronize()
#endif
   kin_local=kin
   deallocate(points,kin)
end subroutine turbine_point_kinematics_gpu
end module m_turbine_point_kinematics_gpu
