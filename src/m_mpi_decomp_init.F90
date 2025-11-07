module m_mpi_decomp_init
#ifdef MPI
   use mpi

   implicit none
   integer :: mpi_rank = 0, mpi_nprocs = 1
   integer :: north = -1, south = -1
   integer :: j_start = 1, j_end = 0
   logical :: periodic_j = .false.
contains
   subroutine mpi_decomp_init(periodic_j_in)
      use mod_dimensions, only : ny,nyg
      implicit none
      logical, intent(in) :: periodic_j_in
      integer :: ierr

      periodic_j = periodic_j_in

      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, mpi_nprocs, ierr)

      if (mod(nyg, mpi_nprocs) /= 0) then
         if (mpi_rank == 0) print *,'Error: ny not divisible by MPI size'
         call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
      end if

      
      j_start = mpi_rank*ny + 1
      j_end   = j_start + ny - 1

      if (periodic_j) then
         north = mod(mpi_rank+1, mpi_nprocs)
         south = mod(mpi_rank-1+mpi_nprocs, mpi_nprocs)
      else
         north = merge(mpi_rank+1, MPI_PROC_NULL, mpi_rank+1 < mpi_nprocs)
         south = merge(mpi_rank-1, MPI_PROC_NULL, mpi_rank-1 >= 0)
      end if
      print '(a,i3,a,i3,a,i5,a,i5,a,i5)', 'rank=',mpi_rank,' nproc=',mpi_nprocs, &
          ' ny=',ny,' j_start=',j_start,' j_end=',j_end

#endif
   end subroutine
end module
