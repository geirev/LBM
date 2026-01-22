module m_mpi_decomp_init
#ifdef MPI
   use mpi
#endif
   implicit none

   integer :: mpi_rank   = 0
   integer :: mpi_nprocs = 1
   integer :: north      = -1
   integer :: south      = -1
   integer :: j_start    = 1
   integer :: j_end      = 0
   logical :: periodic_j = .false.

contains

   subroutine mpi_decomp_init(periodic_j_in)
#ifdef _CUDA
      use cudafor
#endif
      use mod_dimensions, only : ny, nyg, ntiles
      implicit none

      logical, intent(in) :: periodic_j_in
      integer :: ierr

#ifdef MPI
      character(len=MPI_MAX_PROCESSOR_NAME) :: host
      integer :: hostlen
      integer :: shmcomm, shm_rank, shm_size
#endif

#ifdef _CUDA
      integer :: ngpu, dev
#endif

      periodic_j = periodic_j_in

#ifdef MPI
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank,   ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpi_nprocs, ierr)

      ! Hostname for diagnostics
      call MPI_Get_processor_name(host, hostlen, ierr)

      ! Per-node "local rank" (robust GPU selection even if global ranks are shuffled)
      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shmcomm, ierr)
      call MPI_Comm_rank(shmcomm, shm_rank, ierr)
      call MPI_Comm_size(shmcomm, shm_size, ierr)

      ! 1) Enforce that MPI layout matches compile-time tiling
      if (mpi_nprocs /= ntiles) then
         if (mpi_rank == 0) then
            write(*,*) 'Error: mpi_nprocs =', mpi_nprocs, &
                       ' must equal ntiles =', ntiles
         end if
         call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
      end if

      ! 2) Extra sanity check
      if (nyg /= ntiles*ny) then
         if (mpi_rank == 0) then
            write(*,*) 'Error: nyg /= ntiles*ny : ', nyg, ntiles*ny
         end if
         call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
      end if

      ! Print rank placement
      write(*,'(a,i0,a,i0,a,i0,a,a)') 'Rank ', mpi_rank, ' local ', shm_rank, &
           ' of ', shm_size, ' on ', trim(host(1:hostlen))
      call flush(6)
#endif

#ifdef MPI
#ifdef _CUDA
      ierr = cudaGetDeviceCount(ngpu)
      if (ierr /= 0 .or. ngpu <= 0) then
         write(*,*) "cudaGetDeviceCount failed or no GPUs, ierr=", ierr
         call flush(6)
         call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if

      ! Use per-node local rank to choose GPU deterministically on each node
      dev  = mod(shm_rank, ngpu)
      ierr = cudaSetDevice(dev)

      write(*,'(a,i0,a,i0,a,i0,a,a)') 'RANK ', mpi_rank, ' using GPU ', dev, &
           ' of ', ngpu, ' on ', trim(host(1:hostlen))
      write(*,*) "Number of GPUs visible per node: ", ngpu
      write(*,*) "Rank", mpi_rank, "using GPU", dev
      call flush(6)

      if (ierr /= 0) then
         write(*,*) "Rank", mpi_rank, "cudaSetDevice(",dev,") failed, ierr=", ierr
         call flush(6)
         call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
#endif
#endif

      ! Local j-range (ny is tile size)
      j_start = mpi_rank*ny + 1
      j_end   = j_start + ny - 1

#ifdef MPI
      if (periodic_j) then
         north = mod(mpi_rank+1,            mpi_nprocs)
         south = mod(mpi_rank-1+mpi_nprocs, mpi_nprocs)
      else
         north = merge(mpi_rank+1, MPI_PROC_NULL, mpi_rank+1 < mpi_nprocs)
         south = merge(mpi_rank-1, MPI_PROC_NULL, mpi_rank-1 >= 0)
      end if

      print '(a,i3,a,i3,a,i5,a,i5,a,i5)', 'rank=',mpi_rank,' nproc=',mpi_nprocs, &
          ' ny=',ny,' j_start=',j_start,' j_end=',j_end
      call flush(6)

      ! Free the shared-memory communicator
      call MPI_Comm_free(shmcomm, ierr)
#endif

   end subroutine mpi_decomp_init

end module m_mpi_decomp_init
