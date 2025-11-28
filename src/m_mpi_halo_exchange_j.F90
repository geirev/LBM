module m_mpi_halo_exchange_j
contains
subroutine mpi_halo_exchange_j(f,nl)
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
   use m_mpi_decomp_init, only : north, south, mpi_nprocs, mpi_rank
   implicit none
   integer, intent(in) :: nl
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

   call mpi_halo_buffers_alloc(nl)
! Everything is handled in the boundarycond routine for serial code
#ifdef MPI
   if (mpi_nprocs == 1) return
#endif
#ifndef MPI
   return
#endif

   !print *,'mpi: ', mpi_rank, south, north, mpi_nprocs, MPI_PROC_NULL

#ifdef _CUDA
   B = dim3(256,1,1)
   G = dim3( ((nx+2)+B%x-1)/B%x, ((nz+2)+B%y-1)/B%y, (nl+B%z-1)/B%z )
#endif

   call mpi_pack_jplane &
#ifdef _CUDA
   &<<<G,B>>>&
#endif
   (f, 1,  snd_south, nl)

   call mpi_pack_jplane &
#ifdef _CUDA
   &<<<G,B>>>&
#endif
   (f, ny, snd_north, nl)

   ! datatype must match your REAL kind
   if (kind(snd_south(1)) == kind(1.0d0)) then
      mpi_rtype = MPI_DOUBLE_PRECISION
   else
      mpi_rtype = MPI_REAL
   end if

   ! count must be default integer
   count_i = int(plane_elems, kind=4)

#ifdef _CUDA
   istat = cudaDeviceSynchronize()
#endif

   if (south /= MPI_PROC_NULL) then
      call MPI_Irecv(rcv_south, count_i, mpi_rtype, south, 201, MPI_COMM_WORLD, req(1), ierr)
      call MPI_Isend(snd_south, count_i, mpi_rtype, south, 200, MPI_COMM_WORLD, req(2), ierr)
   else
      req(1)=MPI_REQUEST_NULL; req(2)=MPI_REQUEST_NULL
   end if

   if (north /= MPI_PROC_NULL) then
      call MPI_Irecv(rcv_north, count_i, mpi_rtype, north, 200, MPI_COMM_WORLD, req(3), ierr)
      call MPI_Isend(snd_north, count_i, mpi_rtype, north, 201, MPI_COMM_WORLD, req(4), ierr)
   else
      req(3)=MPI_REQUEST_NULL; req(4)=MPI_REQUEST_NULL
   end if

   call MPI_Waitall(4, req, MPI_STATUSES_IGNORE, ierr)

#ifdef _CUDA
   istat = cudaDeviceSynchronize()
#endif

!!!!!!!!!!!!!!!!!!! South unpack
   if (south /= MPI_PROC_NULL) then
      call mpi_unpack_jplane &
#ifdef _CUDA
      &<<<G,B>>>&
#endif
      (f, 0, rcv_south, nl)
   end if

!!!!!!!!!!!!!!!!!!! North unpack
   if (north /= MPI_PROC_NULL) then
      call mpi_unpack_jplane &
#ifdef _CUDA
      &<<<G,B>>>&
#endif
      (f, ny+1, rcv_north, nl)
   endif

#ifdef _CUDA
   istat = cudaDeviceSynchronize(); if (istat/=0) print*,'unpack sync err=',istat
#endif
   call mpi_halo_buffers_free()

end subroutine
end module

