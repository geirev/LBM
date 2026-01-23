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

      ! For global overview
      character(len=MPI_MAX_PROCESSOR_NAME), allocatable :: hosts(:)
      integer, allocatable :: node_ids(:)
      integer :: my_node_id, nnodes, i
      integer :: ngpu_node

      ! Gather per-rank info for one-line-per-rank summary
      integer, allocatable :: all_dev(:), all_j0(:), all_j1(:)
      integer :: j0, j1
         integer :: j
         logical :: seen
#endif

#ifdef _CUDA
      integer :: ngpu, dev
#else
      integer :: dev
      dev = -1
#endif

      periodic_j = periodic_j_in

#ifdef MPI
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank,   ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpi_nprocs, ierr)

      call MPI_Get_processor_name(host, hostlen, ierr)

      call MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shmcomm, ierr)
      call MPI_Comm_rank(shmcomm, shm_rank, ierr)
      call MPI_Comm_size(shmcomm, shm_size, ierr)

      ! 1) Enforce that MPI layout matches compile-time tiling
      if (mpi_nprocs /= ntiles) then
         if (mpi_rank == 0) then
            write(*,*) 'Error: mpi_nprocs =', mpi_nprocs, ' must equal ntiles =', ntiles
         end if
         call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
      end if

      ! 2) Sanity check
      if (nyg /= ntiles*ny) then
         if (mpi_rank == 0) then
            write(*,*) 'Error: nyg /= ntiles*ny : ', nyg, ntiles*ny
         end if
         call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
      end if
#endif

#ifdef MPI
#ifdef _CUDA
      ierr = cudaGetDeviceCount(ngpu)
      if (ierr /= 0 .or. ngpu <= 0) then
         write(*,*) "cudaGetDeviceCount failed or no GPUs, ierr=", ierr
         call flush(6)
         call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if

      dev  = mod(shm_rank, ngpu)
      ierr = cudaSetDevice(dev)
      if (ierr /= 0) then
         write(*,*) "Rank", mpi_rank, "cudaSetDevice(",dev,") failed, ierr=", ierr
         call flush(6)
         call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
      end if
#else
      dev = -1
      ngpu = 0
#endif
#endif

      ! Local j-range (ny is tile size)
      j_start = mpi_rank*ny + 1
      j_end   = j_start + ny - 1

#ifdef MPI
      if (periodic_j) then
         north = mod(mpi_rank+1,              mpi_nprocs)
         south = mod(mpi_rank-1+mpi_nprocs,   mpi_nprocs)
      else
         north = merge(mpi_rank+1, MPI_PROC_NULL, mpi_rank+1 < mpi_nprocs)
         south = merge(mpi_rank-1, MPI_PROC_NULL, mpi_rank-1 >= 0)
      end if

      ! -------------------------------
      ! Global overview (rank 0 prints)
      ! -------------------------------
      allocate(hosts(mpi_nprocs))
      call MPI_Allgather(host, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, &
                         hosts, MPI_MAX_PROCESSOR_NAME, MPI_CHARACTER, MPI_COMM_WORLD, ierr)

      allocate(node_ids(mpi_nprocs))
      node_ids = -1
      nnodes   = 0

      ! Compute each rank's node id deterministically from hostname list (on every rank)
      do i=1, mpi_nprocs
         seen = .false.
         do j=1, i-1
            if (trim(hosts(j)) == trim(hosts(i))) then
               node_ids(i) = node_ids(j)
               seen = .true.
               exit
            end if
         end do
         if (.not. seen) then
            node_ids(i) = nnodes
            nnodes = nnodes + 1
         end if
      end do
      my_node_id = node_ids(mpi_rank+1)

      ! GPUs per node: assume homogeneous nodes, get ngpu from local rank 0 on each node,
      ! then max across world (safe, simple).
#ifdef _CUDA
      ngpu_node = merge(ngpu, 0, shm_rank == 0)
#else
      ngpu_node = 0
#endif
      call MPI_Allreduce(ngpu_node, ngpu_node, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)

      if (mpi_rank == 0) then
         write(*,'(5(a,i0))') 'MPI overview: ranks=', mpi_nprocs,' nyg=', nyg, ' ny=', ny, &
              '  nodes=', nnodes, '  gpus_per_node=', ngpu_node
      end if

      ! --------------------------------------
      ! One line per rank (rank 0 prints all)
      ! --------------------------------------
      allocate(all_dev(mpi_nprocs), all_j0(mpi_nprocs), all_j1(mpi_nprocs))
      j0 = j_start; j1 = j_end
      call MPI_Gather(dev, 1, MPI_INTEGER, all_dev, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(j0,  1, MPI_INTEGER, all_j0,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_Gather(j1,  1, MPI_INTEGER, all_j1,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (mpi_rank == 0) then
         write(*,'(a)') 'Rank  Node   GPU   j_start     j_end  Host'
         do i=1, mpi_nprocs
            write(*,'(i4,2x,i4,2x,i4,2x,i8,2x,i8,2x,a)') &
                 i-1, node_ids(i), all_dev(i), all_j0(i), all_j1(i), trim(hosts(i))
         end do
         call flush(6)
      end if

      deallocate(hosts, node_ids, all_dev, all_j0, all_j1)

      call MPI_Comm_free(shmcomm, ierr)
#endif

   end subroutine mpi_decomp_init

end module m_mpi_decomp_init

