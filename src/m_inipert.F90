module m_inipert
contains
subroutine inipert(rho,u,v,w,uvel,ir)
   use mod_dimensions, only : nx, ny, nyg, nz
   use m_readinfile,  only : rho0, udir, ntx, nty, ntz
#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : mpi_rank, mpi_nprocs
#endif
#ifdef _CUDA
   use cudafor
#endif
   implicit none

   ! Arguments
   real,    intent(inout) :: rho(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(inout) :: u(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(inout) :: v(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(inout) :: w(0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in)    :: uvel(nz)
   integer, intent(in)    :: ir

#ifdef _CUDA
   attributes(device) :: rho, u, v, w, uvel
#endif

   ! Parameters
   real, parameter :: stddev = 1.0e-5
   real, parameter :: pi     = 3.1415927410125732

   ! Host working arrays
   real, allocatable :: rho_global(:,:,:)
   real, allocatable :: rho_local(:,:,:)
#ifdef MPI
   real, allocatable :: sendbuf(:), recvbuf(:)
   integer :: ierr, count, idx
#endif
   integer :: i,j,k

   !--------------------------------------------------------------
   ! 1. Compute random perturbations on the global grid (host)
   !--------------------------------------------------------------
#ifdef MPI
   if (mpi_rank == 0) then
      allocate(rho_global(nx, nyg, nz))
      call random_number(rho_global)
   end if
#else
   ! Serial: we expect ir==0 for the first (and only) call
   if (ir == 0) then
      allocate(rho_global(nx, nyg, nz))
      call random_number(rho_global)
   else
      stop 'inipert (serial): ir /= 0 not supported'
   end if
#endif

   allocate(rho_local(nx, ny, nz))

   !--------------------------------------------------------------
   ! 2. Distribute global field along j (MPI) or copy (serial)
   !--------------------------------------------------------------
#ifdef MPI
   count = nx * ny * nz

   if (mpi_rank == 0) then
      allocate(sendbuf(nx*nyg*nz))
      idx = 0
      do k = 1, nz
         do j = 1, nyg
            do i = 1, nx
               idx = idx + 1
               sendbuf(idx) = rho_global(i,j,k)
            end do
         end do
      end do
   end if

   allocate(recvbuf(count))
   call MPI_Scatter(sendbuf, count, MPI_REAL, &
                    recvbuf, count, MPI_REAL, &
                    0, MPI_COMM_WORLD, ierr)

   idx = 0
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            idx = idx + 1
            rho_local(i,j,k) = recvbuf(idx)
         end do
      end do
   end do

   if (mpi_rank == 0) then
      deallocate(sendbuf, rho_global)
   end if
   deallocate(recvbuf)

#else
   ! No MPI: nyg == ny, so just copy
   rho_local(:,:,:) = rho_global(:,1:ny,:)
   if (allocated(rho_global)) deallocate(rho_global)
#endif

   !--------------------------------------------------------------
   ! 3. Copy local host field into rho interior (1..nx,1..ny,1..nz)
   !--------------------------------------------------------------
   rho(1:nx,1:ny,1:nz) = rho_local(1:nx,1:ny,1:nz)
   deallocate(rho_local)

   !--------------------------------------------------------------
   ! 4. Apply velocity and final rho perturbations on device/host
   !    (this is your original pattern)
   !--------------------------------------------------------------
#ifdef _CUDA
!$cuf kernel do(3) 
#else
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i, j, k) &
!$OMP SHARED(rho, u, v, w, rho0, uvel, stddev, udir, pi)
#endif
   do k=1,nz
      do j=1,ny
         do i=1,nx
            u(i,j,k) = uvel(k)*cos(udir*pi/180.0)
            v(i,j,k) = uvel(k)*sin(udir*pi/180.0)
            w(i,j,k) = 0.0
            rho(i,j,k) = rho0 + stddev*rho(i,j,k)
         end do
      end do
   end do

end subroutine
end module

