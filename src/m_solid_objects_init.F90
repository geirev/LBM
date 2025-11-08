module m_solid_objects_init
contains
subroutine solid_objects_init(blanking_local, lsolids, experiment, ir)
   use mod_dimensions, only : nx, ny, nyg, nz
   use m_cylinder
   use m_city
   use m_city2
   use m_dump_elevation
#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : mpi_rank, mpi_nprocs
#endif
   implicit none
   logical, intent(inout) :: lsolids
   logical, intent(inout) :: blanking_local(0:nx+1,0:ny+1,0:nz+1)
   character(len=*), intent(in) :: experiment
   integer, intent(in)    :: ir
#ifdef _CUDA
   attributes(device) :: blanking_local
#endif

   logical, allocatable :: blanking_global(:,:,:)

#ifdef MPI
   integer :: ierr
   logical, allocatable :: sendbuf(:), recvbuf(:)
   integer :: i,j,k,r,idx
   integer :: j0, j1, chunk
   logical, allocatable :: blank_host(:,:,:)
#endif

   ! Global mask always on host
   allocate(blanking_global(0:nx+1,0:nyg+1,0:nz+1))
   blanking_global = .false.

   !-----------------------------------------------------------
   ! Build global geometry on ir==0 (and rank 0 in MPI case)
   !-----------------------------------------------------------
   if (ir == 0) then
      select case(trim(experiment))
      case('city')
         call city(blanking_global)
      case('city2')
         call city2(blanking_global)
      case('cylinder')
         call cylinder(blanking_global)
      case('airfoil')
         stop 'needs fix airfoil routine for gpu'
      end select
   end if

   select case(trim(experiment))
   case('city')
      lsolids=.true.
   case('city2')
      lsolids=.true.
   case('cylinder')
      lsolids=.true.
   case('airfoil')
      stop 'needs fix airfoil routine for gpu'
   end select

#ifndef MPI
   !-------------------- SERIAL CASE -------------------------
   if (ny /= nyg) stop 'check ny and nyg: must be equal without MPI'
   blanking_local = blanking_global
   if (lsolids) call dump_elevation(blanking_local)
   deallocate(blanking_global)
   return
#else
   !-------------------- MPI CASE ----------------------------

   ! Each rank gets nx*ny*nz interior points
   chunk = nx * ny * nz
   allocate(recvbuf(chunk))

   if (mpi_rank == 0) then
      ! pack all tiles per-rank: [ r=0 block ][ r=1 block ] ...
      allocate(sendbuf(chunk * mpi_nprocs))

      do r = 0, mpi_nprocs-1
         j0  = r*ny + 1
         j1  = j0 + ny - 1
         idx = r*chunk

         do k = 1, nz
            do j = j0, j1
               do i = 1, nx
                  idx = idx + 1
                  sendbuf(idx) = blanking_global(i,j,k)
               end do
            end do
         end do
      end do
   end if

   ! Scatter per-tile chunks
   call MPI_Scatter(sendbuf, chunk, MPI_LOGICAL, &
                    recvbuf, chunk, MPI_LOGICAL, &
                    0, MPI_COMM_WORLD, ierr)

   !---------------------------------------------------------------
   ! Unpack received subdomain into local host mirror
   !---------------------------------------------------------------
   ! Use a host mirror to avoid element-wise host->device copies
   allocate(blank_host(0:nx+1,0:ny+1,0:nz+1))
   blank_host = .false.

   idx = 0
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            idx = idx + 1
            blank_host(i,j,k) = recvbuf(idx)
         end do
      end do
   end do

   ! Halos on host
   blank_host(:,0,:)    = .false.
   blank_host(:,ny+1,:) = .false.
   blank_host(0,:,:)    = .false.
   blank_host(nx+1,:,:) = .false.
   blank_host(:,:,0)    = .false.
   blank_host(:,:,nz+1) = .false.

   blanking_local = blank_host

   ! Clean up
   if (mpi_rank == 0) then
      deallocate(blanking_global, sendbuf)
   else
      deallocate(blanking_global)
   end if
   deallocate(recvbuf, blank_host)

#endif  ! MPI

end subroutine
end module

