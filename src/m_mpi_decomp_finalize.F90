module m_mpi_decomp_finalize
#ifdef MPI
contains
   subroutine mpi_decomp_finalize()
      integer :: ierr
      call MPI_Finalize(ierr)
   end subroutine
#endif
end module


