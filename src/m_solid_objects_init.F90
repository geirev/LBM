module m_solid_objects_init
contains

subroutine solid_objects_init(blanking_local, lsolids, experiment, ir)

   use mod_dimensions, only : nx, ny, nyg, nz
   use m_cylinder
   use m_city
   use m_city2
   use m_city3
   use m_read_bathymetry

#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : mpi_rank, mpi_nprocs
#endif

   implicit none

   logical, intent(inout) :: lsolids
   logical, intent(inout) :: blanking_local(0:nx+1,0:ny+1,0:nz+1)
   character(len=*), intent(in) :: experiment
   integer, intent(in) :: ir

#ifdef _CUDA
   attributes(device) :: blanking_local
#endif

   logical, allocatable :: blanking_global(:,:,:)
   logical, allocatable :: blank_host(:,:,:)

   integer :: i,j,k
   integer(8) :: chunk8

#ifdef MPI
   logical, allocatable :: tilebuf(:)
   integer :: ierr
   integer :: r,idx
   integer :: j0,j1
   integer :: chunk
   integer, parameter :: tag_blank = 1001
#endif


!-----------------------------------------------------------------------
! Number of interior points in one local MPI tile.
!
! With the present dimensions
!
!    nx = 2752
!    ny = 376
!    nz = 80
!
! this is approximately 82.8 million points.
!-----------------------------------------------------------------------
   chunk8 = int(nx,8)*int(ny,8)*int(nz,8)


#ifdef MPI

!-----------------------------------------------------------------------
! Check assumptions of the present one-dimensional j decomposition.
!
! All MPI ranks are assumed to contain exactly ny points in j, so
!
!       nyg = mpi_nprocs*ny
!
! must hold.
!-----------------------------------------------------------------------
   if (nyg /= mpi_nprocs*ny) then

      if (mpi_rank == 0) then
         write(*,*)
         write(*,*) 'ERROR solid_objects_init: inconsistent MPI decomposition'
         write(*,*) 'nyg              = ',nyg
         write(*,*) 'ny               = ',ny
         write(*,*) 'mpi_nprocs        = ',mpi_nprocs
         write(*,*) 'mpi_nprocs*ny     = ',mpi_nprocs*ny
      endif

      call MPI_Abort(MPI_COMM_WORLD,1,ierr)

   endif


!-----------------------------------------------------------------------
! MPI message counts use default INTEGER in this implementation.
! Make sure the local tile fits in that range.
!-----------------------------------------------------------------------
   if (chunk8 > int(huge(chunk),8)) then

      if (mpi_rank == 0) then
         write(*,*) 'ERROR solid_objects_init: MPI tile too large'
         write(*,*) 'Number of tile points = ',chunk8
         write(*,*) 'Maximum MPI count      = ',huge(chunk)
      endif

      call MPI_Abort(MPI_COMM_WORLD,1,ierr)

   endif

   chunk = int(chunk8)


!-----------------------------------------------------------------------
! Only MPI rank 0 constructs and stores the complete global geometry.
!
! Do not use ir here to select the root: the subsequent MPI
! communication explicitly uses mpi_rank=0 as root.
!-----------------------------------------------------------------------
   if (mpi_rank == 0) then

      allocate(blanking_global(0:nx+1,0:nyg+1,0:nz+1),stat=ierr)

      if (ierr /= 0) then
         write(*,*) 'ERROR solid_objects_init: could not allocate blanking_global'
         write(*,*) 'Global dimensions = ',nx+2,nyg+2,nz+2
         call MPI_Abort(MPI_COMM_WORLD,1,ierr)
      endif

      blanking_global = .false.

   endif


!-----------------------------------------------------------------------
! Build global solid-object geometry on MPI rank 0.
!-----------------------------------------------------------------------
   select case(trim(experiment))

   case('city')
      if (mpi_rank == 0) call city(blanking_global)
      lsolids = .true.

   case('city2')
      if (mpi_rank == 0) call city2(blanking_global)
      lsolids = .true.

   case('city3')
      if (mpi_rank == 0) call city3(blanking_global)
      lsolids = .true.

   case('barcelona_2m')
      if (mpi_rank == 0) call read_bathymetry(blanking_global,'2m')
      lsolids = .true.

   case('barcelona_3m')
      if (mpi_rank == 0) call read_bathymetry(blanking_global,'3m')
      lsolids = .true.

   case('barcelona_4m')
      if (mpi_rank == 0) call read_bathymetry(blanking_global,'4m')
      lsolids = .true.

   case('cylinder')
      if (mpi_rank == 0) call cylinder(blanking_global)
      lsolids = .true.

   case('airfoil')
      if (mpi_rank == 0) then
         write(*,*) 'ERROR: airfoil solid-object routine needs GPU/MPI update'
      endif
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)

   case default
      lsolids = .false.

   end select


!-----------------------------------------------------------------------
! Allocate a host copy of the local mask on every rank.
!
! This permits one bulk host->device assignment at the end rather than
! individual element transfers to blanking_local.
!-----------------------------------------------------------------------
   allocate(blank_host(0:nx+1,0:ny+1,0:nz+1),stat=ierr)

   if (ierr /= 0) then
      write(*,*) 'Rank ',mpi_rank, &
                 ': could not allocate blank_host in solid_objects_init'
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif

   blank_host = .false.


!-----------------------------------------------------------------------
! Temporary packed tile.
!
! Only one local-tile-sized buffer is required on each rank.
!
! Rank 0 packs one global j-slab at a time and either
!
!    - unpacks its own tile locally, or
!    - sends the packed tile to the corresponding MPI rank.
!
! This avoids the original global sendbuf(nx*ny*nz*mpi_nprocs),
! which duplicated essentially the complete global mask on rank 0.
!-----------------------------------------------------------------------
   allocate(tilebuf(chunk),stat=ierr)

   if (ierr /= 0) then
      write(*,*) 'Rank ',mpi_rank, &
                 ': could not allocate tilebuf in solid_objects_init'
      call MPI_Abort(MPI_COMM_WORLD,1,ierr)
   endif


   if (mpi_rank == 0) then

      do r = 0,mpi_nprocs-1

         j0 = r*ny + 1
         j1 = j0 + ny - 1

         idx = 0

         do k = 1,nz
            do j = j0,j1
               do i = 1,nx

                  idx = idx + 1
                  tilebuf(idx) = blanking_global(i,j,k)

               enddo
            enddo
         enddo


!-----------------------------------------------------------------------
! Rank 0 keeps its own tile; all other tiles are sent to their owner.
!-----------------------------------------------------------------------
         if (r == 0) then

            idx = 0

            do k = 1,nz
               do j = 1,ny
                  do i = 1,nx

                     idx = idx + 1
                     blank_host(i,j,k) = tilebuf(idx)

                  enddo
               enddo
            enddo

         else

            call MPI_Send(tilebuf,chunk,MPI_LOGICAL,r,tag_blank, &
                          MPI_COMM_WORLD,ierr)

            if (ierr /= MPI_SUCCESS) then
               write(*,*) 'ERROR solid_objects_init: MPI_Send failed'
               write(*,*) 'Destination rank = ',r
               call MPI_Abort(MPI_COMM_WORLD,ierr,ierr)
            endif

         endif

      enddo

   else


!-----------------------------------------------------------------------
! Receive this rank's packed j-slab.
!-----------------------------------------------------------------------
      call MPI_Recv(tilebuf,chunk,MPI_LOGICAL,0,tag_blank, &
                    MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

      if (ierr /= MPI_SUCCESS) then
         write(*,*) 'Rank ',mpi_rank, &
                    ': MPI_Recv failed in solid_objects_init'
         call MPI_Abort(MPI_COMM_WORLD,ierr,ierr)
      endif


!-----------------------------------------------------------------------
! Unpack received interior points.
!-----------------------------------------------------------------------
      idx = 0

      do k = 1,nz
         do j = 1,ny
            do i = 1,nx

               idx = idx + 1
               blank_host(i,j,k) = tilebuf(idx)

            enddo
         enddo
      enddo

   endif


!-----------------------------------------------------------------------
! Local ghost cells do not contain solid objects.
!-----------------------------------------------------------------------
   blank_host(:,0,:)    = .false.
   blank_host(:,ny+1,:) = .false.

   blank_host(0,:,:)    = .false.
   blank_host(nx+1,:,:) = .false.

   blank_host(:,:,0)    = .false.
   blank_host(:,:,nz+1) = .false.


!-----------------------------------------------------------------------
! One bulk host -> device copy in CUDA builds.
! In CPU builds this is a normal array assignment.
!-----------------------------------------------------------------------
   blanking_local = blank_host


!-----------------------------------------------------------------------
! Clean up host memory.
!-----------------------------------------------------------------------
   deallocate(tilebuf)
   deallocate(blank_host)

   if (mpi_rank == 0) then
      deallocate(blanking_global)
   endif


#else

!=======================================================================
! SERIAL CASE
!=======================================================================

!-----------------------------------------------------------------------
! In a serial run ny and nyg are expected to describe the same domain.
!-----------------------------------------------------------------------
   if (ny /= nyg) then
      write(*,*) 'ERROR solid_objects_init: serial run requires ny = nyg'
      write(*,*) 'ny  = ',ny
      write(*,*) 'nyg = ',nyg
      error stop
   endif


!-----------------------------------------------------------------------
! Allocate complete solid mask.
!-----------------------------------------------------------------------
   allocate(blanking_global(0:nx+1,0:nyg+1,0:nz+1))

   blanking_global = .false.


!-----------------------------------------------------------------------
! Build geometry.
!-----------------------------------------------------------------------
   select case(trim(experiment))

   case('city')
      call city(blanking_global)
      lsolids = .true.

   case('city2')
      call city2(blanking_global)
      lsolids = .true.

   case('city3')
      call city3(blanking_global)
      lsolids = .true.

   case('barcelona_2m')
      call read_bathymetry(blanking_global,'2m')
      lsolids = .true.

   case('barcelona_3m')
      call read_bathymetry(blanking_global,'3m')
      lsolids = .true.

   case('barcelona_4m')
      call read_bathymetry(blanking_global,'4m')
      lsolids = .true.

   case('cylinder')
      call cylinder(blanking_global)
      lsolids = .true.

   case('airfoil')
      error stop 'needs fix airfoil routine for gpu'

   case default
      lsolids = .false.

   end select


!-----------------------------------------------------------------------
! Copy global mask into local mask.
!-----------------------------------------------------------------------
   blanking_local = blanking_global

   deallocate(blanking_global)

#endif


end subroutine solid_objects_init

end module m_solid_objects_init
