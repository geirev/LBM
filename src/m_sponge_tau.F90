module m_sponge_tau
contains
subroutine sponge_tau(tau,tau_max,nsponge_i,nsponge_j)
   ! Ramps the current (turbulence-scheme-computed) tau field up
   ! toward tau_max near open boundaries. Call this every timestep,
   ! after your turbulence scheme updates tau and before the
   ! collision step - it modifies tau in place.
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile,  only : ibnd,jbnd
#ifdef _CUDA
   use m_readinfile,  only : ntx,nty,ntz
#endif
#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : north, south
#endif
   use m_wtime
   use m_sponge_tau_kernel
   implicit none
   real, intent(inout) :: tau(0:nx+1,0:ny+1,0:nz+1)
   real, value           :: tau_max
   integer, value        :: nsponge_i
   integer, value        :: nsponge_j
#ifdef _CUDA
   attributes(device) :: tau
   integer :: tx, ty, tz, bx, by, bz
#endif
   integer, parameter :: icpu=13
   logical :: j0_is_phys, jN_is_phys

   call cpustart()

#ifdef MPI
   j0_is_phys = (south == MPI_PROC_NULL)
   jN_is_phys = (north == MPI_PROC_NULL)
#else
   j0_is_phys = .true.
   jN_is_phys = .true.
#endif

#ifdef _CUDA
   tx = ntx; bx = (nx + tx - 1) / tx
   ty = nty; by = (ny + ty - 1) / ty
   tz = ntz; bz = (nz + tz - 1) / tz
#endif
   call sponge_tau_kernel &
#ifdef _CUDA
        <<<dim3(bx,by,bz), dim3(tx,ty,tz)>>> &
#endif
        (tau,tau_max,nsponge_i,nsponge_j, &
         ibnd,jbnd,j0_is_phys,jN_is_phys)

   call cpufinish(icpu)
end subroutine
end module
