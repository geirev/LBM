module m_dump_elevation
contains
subroutine dump_elevation(lblanking)
   use mod_dimensions
   use m_tecfld
   use m_readinfile, only : lmeasurements
   implicit none
   logical, intent(in)  :: lblanking(0:nx+1,0:ny+1,0:nz+1)
   logical, allocatable :: lblanking_h(:,:,:)
   real, allocatable    :: elevation(:,:)
   integer(kind=4), allocatable    :: elevation3(:,:,:)
   integer i,j,k,iunit
#ifdef _CUDA
   attributes(device) :: lblanking
#endif

   allocate(lblanking_h(0:nx+1,0:ny+1,0:nz+1))
   lblanking_h=lblanking

   allocate(elevation(nx,ny))
   allocate(elevation3(nx,ny,nz))
   elevation=0.0
   elevation3=0

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
   call tecfld('elevation',nx,ny,1,elevation)
   if (lmeasurements) then
      open(newunit=iunit,file='blanking3D.uf',form="unformatted", status='unknown')
         write(iunit)elevation3
      close(iunit)
   endif

   deallocate(lblanking_h,elevation,elevation3)
end subroutine
end module
