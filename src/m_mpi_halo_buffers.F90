module m_mpi_halo_buffers
#ifdef MPI
   use mpi
#endif
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions, only: nx, ny, nz
   use m_mpi_decomp_init,  only: north, south
   implicit none
   ! Contiguous device buffers for one full j-plane
   real, allocatable :: snd_south(:), rcv_south(:), snd_north(:), rcv_north(:)
#ifdef _CUDA
   attributes(device)  :: snd_south
   attributes(device)  :: rcv_south
   attributes(device)  :: snd_north
   attributes(device)  :: rcv_north
#endif
 !  integer(kind=8) :: plane_elems = 0
contains
   subroutine mpi_halo_buffers_alloc(nl)
      use mod_dimensions, only: nx,nz
      implicit none
      integer, intent(in) :: nl
      integer :: plane_elems
      plane_elems = nl*(nx+2)*(nz+2)
      if (.not.allocated(snd_south)) allocate(snd_south(plane_elems))
      if (.not.allocated(rcv_south)) allocate(rcv_south(plane_elems))
      if (.not.allocated(snd_north)) allocate(snd_north(plane_elems))
      if (.not.allocated(rcv_north)) allocate(rcv_north(plane_elems))
   end subroutine
   subroutine mpi_halo_buffers_free()
      implicit none
      if (allocated(snd_south)) deallocate(snd_south)
      if (allocated(rcv_south)) deallocate(rcv_south)
      if (allocated(snd_north)) deallocate(snd_north)
      if (allocated(rcv_north)) deallocate(rcv_north)
   end subroutine
end module

