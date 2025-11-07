module m_mpi_halo_exchange_j
contains
subroutine mpi_halo_exchange_j(f)
#ifdef MPI
   use mpi
   use m_mpi_pack_jplane
   use m_mpi_unpack_jplane
   use m_mpi_halo_buffers           ! snd_*/rcv_* and plane_elems
#endif
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only : nx, ny, nz
   use mod_D3Q27setup, only : nl
   use m_mpi_decomp_init, only : north, south
   implicit none
   real :: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: f
   type(dim3) :: B, G
   integer :: istat
#endif
#ifdef MPI
   integer :: ierr, req(4)
   integer :: count_i, mpi_rtype
#endif

#ifdef _CUDA
   B = dim3(16,8,4)
   G = dim3( ((nx+2)+B%x-1)/B%x, ((nz+2)+B%y-1)/B%y, (nl+B%z-1)/B%z )
#endif

   call mpi_pack_jplane &
#ifdef _CUDA
   &<<<G,B>>>&
#endif
   (f, 1,  snd_south)

   call mpi_pack_jplane &
#ifdef _CUDA
   &<<<G,B>>>&
#endif
   (f, ny, snd_north)

#ifdef MPI
   ! datatype must match your REAL kind
   if (kind(snd_south(1)) == kind(1.0d0)) then
      mpi_rtype = MPI_DOUBLE_PRECISION
   else
      mpi_rtype = MPI_REAL
   end if

   ! count must be default integer
   count_i = int(plane_elems, kind=4)

   ! recv what we need for our ghosts
   call MPI_Irecv(rcv_south, count_i, mpi_rtype, south, 201, MPI_COMM_WORLD, req(1), ierr)
   call MPI_Irecv(rcv_north, count_i, mpi_rtype, north, 200, MPI_COMM_WORLD, req(3), ierr)

   ! send our boundaries out
   call MPI_Isend(snd_south, count_i, mpi_rtype, south, 200, MPI_COMM_WORLD, req(2), ierr)
   call MPI_Isend(snd_north, count_i, mpi_rtype, north, 201, MPI_COMM_WORLD, req(4), ierr)

   call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, ierr)
#endif

   call mpi_unpack_jplane &
#ifdef _CUDA
   &<<<G,B>>>&
#endif
   (f, 0,    rcv_south)

   call mpi_unpack_jplane &
#ifdef _CUDA
   &<<<G,B>>>&
#endif
   (f, ny+1, rcv_north)

#ifdef _CUDA
   istat = cudaDeviceSynchronize(); if (istat/=0) print*,'unpack sync err=',istat
#endif

end subroutine
end module

