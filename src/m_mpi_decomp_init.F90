module m_mpi_decomp_init
#ifdef MPI
   use mpi
   implicit none
   integer :: mpi_rank = 0, mpi_nprocs = 1
   integer :: north = MPI_PROC_NULL, south = MPI_PROC_NULL
   integer :: ny_local = 0, j_start = 1, j_end = 0
   logical :: periodic_j = .false.
contains
   subroutine mpi_decomp_init(ny_global, periodic_j_in)
      integer, intent(in) :: ny_global
      logical, intent(in) :: periodic_j_in
      integer :: ierr

      periodic_j = periodic_j_in

      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpi_nprocs, ierr)

      if (mod(ny_global, mpi_nprocs) /= 0) then
         if (mpi_rank == 0) print *,'Error: ny not divisible by MPI size'
         call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
      end if

      ny_local = ny_global / mpi_nprocs
      j_start  = mpi_rank*ny_local + 1
      j_end    = j_start + ny_local - 1

      if (periodic_j) then
         north = mod(mpi_rank+1, mpi_nprocs)
         south = mod(mpi_rank-1+mpi_nprocs, mpi_nprocs)
      else
         north = merge(mpi_rank+1, MPI_PROC_NULL, mpi_rank+1 < mpi_nprocs)
         south = merge(mpi_rank-1, MPI_PROC_NULL, mpi_rank-1 >= 0)
      end if
      print '(7(a,i3))','mpi_rank=',mpi_rank,' mpi_nprocs=',mpi_nprocs,' ny_local=',ny_local,&
              ' j_start=',j_start,' j_end=', j_end,'north=',north,' south=',south
   end subroutine
#endif
end module
