module m_inflow_turbulence_update
contains

subroutine inflow_turbulence_update(uu,vv,ww,rr,nrturb,lfirst)

   use mod_dimensions, only : nx,ny,nz,nyg
   use m_inflow_turbulence_init, only : iturb_pos,nturb_perim
   use m_inflow_turbulence_compute
#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : mpi_rank,mpi_nprocs
#endif

   implicit none

   integer, intent(in) :: nrturb
   logical, intent(in) :: lfirst

   real, intent(inout) :: uu(nturb_perim,nz,0:nrturb)
   real, intent(inout) :: vv(nturb_perim,nz,0:nrturb)
   real, intent(inout) :: ww(nturb_perim,nz,0:nrturb)
   real, intent(inout) :: rr(nturb_perim,nz,0:nrturb)

#ifdef _CUDA
   attributes(device) :: uu,vv,ww,rr
#endif

   integer :: i0,i1,j0g,j1g,ns
   integer :: p,k,it,s

   real, allocatable :: uu_g(:,:,:),vv_g(:,:,:),ww_g(:,:,:),rr_g(:,:,:)

#ifdef MPI
   integer :: ierr,r,jg,joff,np,idx,total
   integer, allocatable :: counts(:),displs(:)
   real, allocatable :: ubuf(:),vbuf(:),wbuf(:),rbuf(:)
#endif

   ! Forcing perimeter: iturb_pos=0 is on the boundary fluid nodes.
   i0  = 1+iturb_pos
   i1  = nx-iturb_pos
   j0g = 1+iturb_pos
   j1g = nyg-iturb_pos

   ns = 2*(i1-i0+1) + 2*(j1g-j0g+1) - 4

#ifndef MPI

   allocate(uu_g(ns,nz,0:nrturb),vv_g(ns,nz,0:nrturb))
   allocate(ww_g(ns,nz,0:nrturb),rr_g(ns,nz,0:nrturb))

   call inflow_turbulence_compute( &
      uu_g,vv_g,ww_g,rr_g,ns,nz,nrturb,lfirst)

   do it=0,nrturb
   do k=1,nz
   do p=1,nturb_perim
      uu(p,k,it)=uu_g(p,k,it)
      vv(p,k,it)=vv_g(p,k,it)
      ww(p,k,it)=ww_g(p,k,it)
      rr(p,k,it)=rr_g(p,k,it)
   enddo
   enddo
   enddo

   deallocate(uu_g,vv_g,ww_g,rr_g)

#else

   allocate(counts(0:mpi_nprocs-1),displs(0:mpi_nprocs-1))

   ! Number of perimeter points owned by each j tile.  Local perimeter
   ! arrays are ordered by increasing global perimeter index s.
   do r=0,mpi_nprocs-1
      joff=r*ny
      np=0
      do s=1,ns
         call perimeter_s_to_j(s,i0,i1,j0g,j1g,jg)
         if (jg > joff .and. jg <= joff+ny) np=np+1
      enddo
      counts(r)=np
   enddo

   ! -----------------------------------------------------------------
   ! Continue the AR(1) sequence across generated blocks.
   ! On all calls after the first, gather the last local realization
   ! from every rank and reconstruct the global perimeter realization
   ! in uu_g(:,:,nrturb), etc. on rank 0.  inflow_turbulence_compute
   ! then copies this nrturb slice to slice 0 before generating the new
   ! innovations 1:nrturb.
   ! -----------------------------------------------------------------
   if (mpi_rank == 0) then
      allocate(uu_g(ns,nz,0:nrturb),vv_g(ns,nz,0:nrturb))
      allocate(ww_g(ns,nz,0:nrturb),rr_g(ns,nz,0:nrturb))
   endif

   if (.not.lfirst) then

      ! Gather one complete (local perimeter,nz) realization per rank.
      ! counts/displs are temporarily expressed in REAL values.
      do r=0,mpi_nprocs-1
         counts(r)=counts(r)*nz
      enddo

      displs(0)=0
      do r=1,mpi_nprocs-1
         displs(r)=displs(r-1)+counts(r-1)
      enddo
      total=displs(mpi_nprocs-1)+counts(mpi_nprocs-1)

      if (mpi_rank == 0) then
         allocate(ubuf(total),vbuf(total),wbuf(total),rbuf(total))
      else
         allocate(ubuf(1),vbuf(1),wbuf(1),rbuf(1))
      endif

      call MPI_Gatherv(uu(1,1,nrturb),counts(mpi_rank),MPI_REAL, &
                       ubuf,counts,displs,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Gatherv(vv(1,1,nrturb),counts(mpi_rank),MPI_REAL, &
                       vbuf,counts,displs,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Gatherv(ww(1,1,nrturb),counts(mpi_rank),MPI_REAL, &
                       wbuf,counts,displs,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_Gatherv(rr(1,1,nrturb),counts(mpi_rank),MPI_REAL, &
                       rbuf,counts,displs,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      if (mpi_rank == 0) then
         ! Undo the rank-wise packing.  Within each rank the local
         ! storage order is p-fastest, then k, matching the loops below.
         do r=0,mpi_nprocs-1
            joff=r*ny
            idx=displs(r)
            do k=1,nz
            do s=1,ns
               call perimeter_s_to_j(s,i0,i1,j0g,j1g,jg)
               if (jg > joff .and. jg <= joff+ny) then
                  idx=idx+1
                  uu_g(s,k,nrturb)=ubuf(idx)
                  vv_g(s,k,nrturb)=vbuf(idx)
                  ww_g(s,k,nrturb)=wbuf(idx)
                  rr_g(s,k,nrturb)=rbuf(idx)
               endif
            enddo
            enddo
         enddo
      endif

      deallocate(ubuf,vbuf,wbuf,rbuf)

      ! Restore counts to number of perimeter points for use below.
      do r=0,mpi_nprocs-1
         counts(r)=counts(r)/nz
      enddo

   endif

   if (mpi_rank == 0) then
      call inflow_turbulence_compute( &
         uu_g,vv_g,ww_g,rr_g,ns,nz,nrturb,lfirst)
   endif

   ! -----------------------------------------------------------------
   ! Scatter the newly generated complete time block.  Each rank gets
   ! only the perimeter points lying in its local j tile.
   ! -----------------------------------------------------------------
   do r=0,mpi_nprocs-1
      counts(r)=counts(r)*nz*(nrturb+1)
   enddo

   displs(0)=0
   do r=1,mpi_nprocs-1
      displs(r)=displs(r-1)+counts(r-1)
   enddo
   total=displs(mpi_nprocs-1)+counts(mpi_nprocs-1)

   if (mpi_rank == 0) then
      allocate(ubuf(total),vbuf(total),wbuf(total),rbuf(total))

      do r=0,mpi_nprocs-1
         joff=r*ny
         idx=displs(r)
         do it=0,nrturb
         do k=1,nz
         do s=1,ns
            call perimeter_s_to_j(s,i0,i1,j0g,j1g,jg)
            if (jg > joff .and. jg <= joff+ny) then
               idx=idx+1
               ubuf(idx)=uu_g(s,k,it)
               vbuf(idx)=vv_g(s,k,it)
               wbuf(idx)=ww_g(s,k,it)
               rbuf(idx)=rr_g(s,k,it)
            endif
         enddo
         enddo
         enddo
      enddo
   else
      allocate(ubuf(1),vbuf(1),wbuf(1),rbuf(1))
   endif

   call MPI_Scatterv(ubuf,counts,displs,MPI_REAL, &
                     uu,counts(mpi_rank),MPI_REAL,0,MPI_COMM_WORLD,ierr)
   call MPI_Scatterv(vbuf,counts,displs,MPI_REAL, &
                     vv,counts(mpi_rank),MPI_REAL,0,MPI_COMM_WORLD,ierr)
   call MPI_Scatterv(wbuf,counts,displs,MPI_REAL, &
                     ww,counts(mpi_rank),MPI_REAL,0,MPI_COMM_WORLD,ierr)
   call MPI_Scatterv(rbuf,counts,displs,MPI_REAL, &
                     rr,counts(mpi_rank),MPI_REAL,0,MPI_COMM_WORLD,ierr)

   deallocate(ubuf,vbuf,wbuf,rbuf,counts,displs)

   if (mpi_rank == 0) then
      deallocate(uu_g,vv_g,ww_g,rr_g)
   endif

#endif

contains

subroutine perimeter_s_to_j(s,i0,i1,j0,j1,j)

   implicit none

   integer, intent(in)  :: s,i0,i1,j0,j1
   integer, intent(out) :: j

   integer :: nxp,nyp,q

   nxp=i1-i0+1
   nyp=j1-j0+1

   if (s <= nxp) then
      j=j0
   elseif (s <= nxp+nyp-1) then
      q=s-nxp
      j=j0+q
   elseif (s <= nxp+nyp-1+nxp-1) then
      j=j1
   else
      q=s-(nxp+nyp-1+nxp-1)
      j=j1-q
   endif

end subroutine perimeter_s_to_j

end subroutine inflow_turbulence_update

end module m_inflow_turbulence_update
