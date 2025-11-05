module m_mpi_halo_j
#ifdef MPI
   use mpi
#endif
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only: nx, nz, nl
   use m_mpi_decomp_init,  only: ny_local, north, south
   implicit none
   ! Contiguous device buffers for one full j-plane
   real, allocatable :: snd_south(:), rcv_south(:), snd_north(:), rcv_north(:)
#ifdef _CUDA
   attributes(device)  :: snd_south
   attributes(device)  :: rcv_south
   attributes(device)  :: snd_north
   attributes(device)  :: rcv_north
#endif
   integer(kind=8) :: plane_elems = 0
contains
   subroutine halo_buffers_alloc()
      plane_elems = int(nl,8)*int(nx+2,8)*int(nz+2,8)
      allocate(snd_south(plane_elems), rcv_south(plane_elems))
      allocate(snd_north(plane_elems), rcv_north(plane_elems))
   end subroutine
   subroutine halo_buffers_free()
      if (allocated(snd_south)) deallocate(snd_south, rcv_south, snd_north, rcv_north)
   end subroutine

#ifdef _CUDA
   attributes(global)&
#endif
   subroutine pack_jplane(f, j_plane, buf)
      integer, value :: j_plane
      real :: f(:,:,:,:)
      real :: buf(:)
#ifdef _CUDA
      attributes(device)  :: f
      attributes(device)  :: buf
#endif
      integer :: i,k,l
      integer(kind=8) :: idx, nxi
#ifdef _CUDA
      i = (blockIdx%x-1)*blockDim%x + threadIdx%x
      k = (blockIdx%y-1)*blockDim%y + threadIdx%y
      l = (blockIdx%z-1)*blockDim%z + threadIdx%z
      if (i<0 .or. i>nx+1) return
      if (k<0 .or. k>nz+1) return
      if (l<1 .or. l>nl  ) return
#else
      do k=0,nz+1
      do i=0,nx+1
      do l=1,nl
#endif
         nxi = int(nx+2,8)
         idx = ( int(k,8)*nxi + int(i,8) )*int(nl,8) + int(l,8)
         buf(idx) = f(l,i,j_plane,k)
#ifndef _CUDA
      enddo
      enddo
      enddo
#endif
   end subroutine

#ifdef _CUDA
   attributes(global)&
#endif
   subroutine unpack_jplane(f, j_plane, buf)
      integer, value :: j_plane
      real :: f(:,:,:,:)
      real :: buf(:)
#ifdef _CUDA
      attributes(device)  :: f
      attributes(device)  :: buf
#endif
      integer :: i,k,l
      integer(kind=8) :: idx, nxi
#ifdef _CUDA
      i = (blockIdx%x-1)*blockDim%x + threadIdx%x
      k = (blockIdx%y-1)*blockDim%y + threadIdx%y
      l = (blockIdx%z-1)*blockDim%z + threadIdx%z
      if (i<0 .or. i>nx+1) return
      if (k<0 .or. k>nz+1) return
      if (l<1 .or. l>nl  ) return
#else
      do k=0,nz+1
      do i=0,nx+1
      do l=1,nl
#endif
      nxi = int(nx+2,8)
      idx = ( int(k,8)*nxi + int(i,8) )*int(nl,8) + int(l,8)
      f(l,i,j_plane,k) = buf(idx)
#ifndef _CUDA
      enddo
      enddo
      enddo
#endif
   end subroutine

   ! CUDA-aware halo exchange for j; call this immediately before DRIFT
   subroutine halo_exchange_j(f)
      real :: f(:,:,:,:)
#ifdef _CUDA
      attributes(device)  :: f
#endif
      integer :: ierr, req(4), stat(MPI_STATUS_SIZE,4)
#ifdef _CUDA
      type(dim3) :: B, G
#endif
      external MPI_Irecv, MPI_Isend, MPI_Irecv, MPI_Isend
      integer istat

#ifdef _CUDA
      B = dim3(16,8,4)
      G = dim3( ((nx+2)+B%x-1)/B%x, ((nz+2)+B%y-1)/B%y, (nl+B%z-1)/B%z )
#endif

      ! Pack j=1 and j=ny_local planes
      call pack_jplane&
#ifdef _CUDA
      &<<<G,B>>>&
#endif
      &(f, 1,          snd_south)

      call pack_jplane&
#ifdef _CUDA
      &<<<G,B>>>&
#endif
      &(f, ny_local,   snd_north)

#ifdef _CUDA
      istat=cudaDeviceSynchronize()
#endif

      ! Nonblocking CUDA-aware MPI: south uses tags 101/102, north mirrored
      call MPI_Irecv(rcv_south, plane_elems, MPI_REAL, south, 102, MPI_COMM_WORLD, req(1), ierr)
      call MPI_Isend(snd_south, plane_elems, MPI_REAL, south, 101, MPI_COMM_WORLD, req(2), ierr)
      call MPI_Irecv(rcv_north, plane_elems, MPI_REAL, north, 101, MPI_COMM_WORLD, req(3), ierr)
      call MPI_Isend(snd_north, plane_elems, MPI_REAL, north, 102, MPI_COMM_WORLD, req(4), ierr)

      call MPI_Waitall(4, req, stat, ierr)

      ! Unpack into ghosts j=0 and j=ny_local+1
      call unpack_jplane&
#ifdef _CUDA
      &<<<G,B>>>&
#endif
      &(f, 0,            rcv_south)

      call unpack_jplane&
#ifdef _CUDA
      &<<<G,B>>>&
#endif
      &(f, ny_local+1,   rcv_north)

#ifdef _CUDA
      istat=cudaDeviceSynchronize()
#endif

   end subroutine
end module

