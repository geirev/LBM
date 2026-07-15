module m_dump_elevation
contains
subroutine dump_elevation(lblanking)
   use mod_dimensions, only : nx,ny,nz
   use m_tecfld
   use m_readinfile, only : lmeasurements
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank,mpi_nprocs
#endif
   implicit none
   logical            :: lblanking(0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: lblanking
#endif

   character(len=4)                :: ctile
   character(len=100)              :: fname
   character(len=100)              :: directory
   logical,         allocatable    :: lblanking_h(:,:,:)
   real,            allocatable    :: elevation(:,:)
   integer(kind=4), allocatable    :: elevation3(:,:,:)
   integer i,j,k,iunit,ir,np

   ir=0
   np=1
#ifdef MPI
   ir=mpi_rank
   np=mpi_nprocs
#endif

   allocate(elevation(nx,ny+1))
   allocate(lblanking_h(0:nx+1,0:ny+1,0:nz+1))
   allocate(elevation3(nx,ny,nz))
   elevation=0.0
   elevation3=0
   lblanking_h=lblanking

   do j=1,ny
   do i=1,nx
   do k=1,nz
      if (lblanking_h(i,j,k)) then
          elevation(i,j)=k
          elevation3(i,j,k)=k
      else
          exit
      endif
   enddo
   enddo
   enddo
   elevation(:,ny+1)=elevation(:,ny)

   directory='output/'
   call system('mkdir -p '//trim(directory))

   write(ctile,'(i4.4)')ir

   fname=trim(directory)//'tec_elevation'//trim(ctile)
   print *,'Elevation: ',trim(fname)
   call tecfld(trim(fname),'elevation',nx,ny,1,elevation,ir,np)

#ifndef MPI
   if (lmeasurements) then
      open(newunit=iunit,file='blanking3D.uf',form="unformatted", status='unknown')
         write(iunit)elevation3
      close(iunit)
   endif
#endif

   deallocate(lblanking_h,elevation,elevation3)
end subroutine
end module
