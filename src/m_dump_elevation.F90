module m_dump_elevation
contains
subroutine dump_elevation(lblanking)
   use mod_dimensions
   use m_tecfld
   implicit none
   logical, intent(in)  :: lblanking(0:nx,0:ny,0:nz)
   logical, allocatable :: lblanking_h(:,:,:)
   real, allocatable :: elevation(:,:)
   integer i,j,k
#ifdef _CUDA
   attributes(device) :: lblanking
#endif

   allocate(lblanking_h(0:nx,0:ny,0:nz))
   lblanking_h=lblanking

   allocate(elevation(nx,ny))
   elevation=0.0

   do j=1,ny
   do i=1,nx
   do k=1,nz
      if (lblanking_h(i,j,k)) then
          elevation(i,j)=k
      else
          exit
      endif
   enddo
   enddo
   enddo
   call tecfld('elevation',nx,ny,1,elevation)
   deallocate(lblanking_h,elevation)
end subroutine
end module
