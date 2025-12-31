module m_turbine_deposit_gpu
contains
subroutine turbine_deposit_gpu(F_turb, xg, yg, zg, Fvec, np)
  use mod_dimensions, only : nx, ny, nz, nyg
#ifdef _CUDA
  use cudafor
#endif
  use m_turbine_deposit_gpu_kernel
#ifdef MPI
  use m_mpi_decomp_init, only : j_start, j_end, mpi_rank
#endif
  implicit none

  integer, intent(in)        :: np
#ifdef _CUDA
  real, device,   intent(inout)     :: F_turb(3,0:nx+1,0:ny+1,0:nz+1)
  real, device,   intent(in)        :: xg(np), yg(np), zg(np), Fvec(3,np)
#else
  real,           intent(inout)     :: F_turb(3,0:nx+1,0:ny+1,0:nz+1)
  real,           intent(in)        :: xg(np), yg(np), zg(np), Fvec(3,np)
#endif

#ifdef _CUDA
  integer, parameter :: tpb = 256   ! threads per block (tune: 128/256/512)
#endif
  integer istat

  real :: epsilon, sigma, sigma2
  integer :: krad

#ifndef _CUDA
  ! ---- CPU fallback (optional): you can keep your original code here ----
  stop "turbine_deposit: compiled without _CUDA, please keep CPU version or add fallback."
#else
  type(dim3) :: grid, block
  integer    :: jstart_in, jend_in

  epsilon = 2.0
  sigma   = epsilon / sqrt(2.0)
  sigma2  = sigma*sigma
  krad    = min(5,int(ceiling(3.0*sigma)))

#ifdef MPI
  jstart_in = j_start
  jend_in   = j_end
#else
  jstart_in = 1
  jend_in   = nyg
#endif

  block = dim3(TPB, 1, 1)
  grid  = dim3(np,  1, 1)

  call turbine_deposit_gpu_kernel<<<grid, block>>>( &
       F_turb, xg, yg, zg, Fvec, np, &
       nx, ny, nz, nyg, jstart_in, jend_in, krad, sigma2)

   istat = cudaDeviceSynchronize()
   if (istat /= cudaSuccess) then
     write(*,*) "Kernel execution failed:", cudaGetErrorString(istat)
     stop
   end if
#endif

end subroutine turbine_deposit_gpu
end module m_turbine_deposit_gpu
