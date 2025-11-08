module m_inflow_turbulence_update
contains
subroutine inflow_turbulence_update(uu,vv,ww,rr,nrturb,lfirst)
   use mod_dimensions, only : ny, nz, nyg
   use m_inflow_turbulence_compute
#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only :  mpi_rank, mpi_nprocs
#endif
   implicit none
   integer, intent(in) :: nrturb
   logical, intent(in) :: lfirst
   real, intent(inout) :: uu(ny,nz,0:nrturb)
   real, intent(inout) :: vv(ny,nz,0:nrturb)
   real, intent(inout) :: ww(ny,nz,0:nrturb)
   real, intent(inout) :: rr(ny,nz,0:nrturb)
#ifdef _CUDA
   attributes(device) :: uu, vv, ww, rr
#endif

   real, allocatable :: uu_g(:,:,:), vv_g(:,:,:), ww_g(:,:,:), rr_g(:,:,:)
#ifdef MPI
   integer :: ierr
   ! packed send buffers on root only
   real, allocatable :: ubuf(:), vbuf(:), wbuf(:), rbuf(:)
   integer :: r, j0, j1, it, k, j, idx, chunk
#endif

#ifndef MPI
   ! ------------------ Serial case: compute directly into local arrays ------------------
      allocate(uu_g(ny,nz,0:nrturb))
      allocate(vv_g(ny,nz,0:nrturb))
      allocate(ww_g(ny,nz,0:nrturb))
      allocate(rr_g(ny,nz,0:nrturb))
      call inflow_turbulence_compute(uu_g,vv_g,ww_g,rr_g,ny,nz,nrturb,lfirst)
      uu=uu_g
      vv=vv_g
      ww=ww_g
      rr=rr_g
      return
#else
   ! ------------------ MPI case: rank 0 generates global, packs, then scatters ---------
   chunk = ny * nz * (nrturb + 1)

   if (mpi_rank == 0) then
      allocate(uu_g(nyg,nz,0:nrturb))
      allocate(vv_g(nyg,nz,0:nrturb))
      allocate(ww_g(nyg,nz,0:nrturb))
      allocate(rr_g(nyg,nz,0:nrturb))

      call inflow_turbulence_compute(uu_g,vv_g,ww_g,rr_g,nyg,nz,nrturb,lfirst)

      allocate(ubuf(chunk*mpi_nprocs), vbuf(chunk*mpi_nprocs), &
               wbuf(chunk*mpi_nprocs), rbuf(chunk*mpi_nprocs))

      ! Pack each rank's (j0:j1, 1:nz, 0:nrturb) tile contiguously: j-fastest, then k, then t
      do r = 0, mpi_nprocs-1
         j0  = r*ny + 1
         j1  = j0 + ny - 1
         idx = r*chunk
         do it = 0, nrturb
            do k = 1, nz
               do j = j0, j1
                  idx = idx + 1
                  ubuf(idx) = uu_g(j,k,it)
               end do
            end do
         end do
      end do

      do r = 0, mpi_nprocs-1
         j0  = r*ny + 1
         j1  = j0 + ny - 1
         idx = r*chunk
         do it = 0, nrturb
            do k = 1, nz
               do j = j0, j1
                  idx = idx + 1
                  vbuf(idx) = vv_g(j,k,it)
               end do
            end do
         end do
      end do

      do r = 0, mpi_nprocs-1
         j0  = r*ny + 1
         j1  = j0 + ny - 1
         idx = r*chunk
         do it = 0, nrturb
            do k = 1, nz
               do j = j0, j1
                  idx = idx + 1
                  wbuf(idx) = ww_g(j,k,it)
               end do
            end do
         end do
      end do

      do r = 0, mpi_nprocs-1
         j0  = r*ny + 1
         j1  = j0 + ny - 1
         idx = r*chunk
         do it = 0, nrturb
            do k = 1, nz
               do j = j0, j1
                  idx = idx + 1
                  rbuf(idx) = rr_g(j,k,it)
               end do
            end do
         end do
      end do
   end if

   ! Scatter equal-size contiguous chunks to each rank’s local 3D arrays
   call MPI_Scatter(ubuf, chunk, MPI_REAL, uu, chunk, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Scatter(vbuf, chunk, MPI_REAL, vv, chunk, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Scatter(wbuf, chunk, MPI_REAL, ww, chunk, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
   call MPI_Scatter(rbuf, chunk, MPI_REAL, rr, chunk, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

   if (mpi_rank == 0) then
      deallocate(uu_g, vv_g, ww_g, rr_g, ubuf, vbuf, wbuf, rbuf)
   end if
#endif

end subroutine inflow_turbulence_update
end module

