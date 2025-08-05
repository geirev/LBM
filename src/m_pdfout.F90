module m_pdfout

contains
subroutine pdfout(filename,variables_string,num_of_variables,lblanking,f)
   use mod_dimensions
   use m_tecplot
   implicit none
   character(len=*), intent(in) :: filename
   character(len=*), intent(in) :: variables_string
   integer, intent(in) :: num_of_variables
   logical, intent(in) :: lblanking(nx,ny,nz)
   real, intent(in)    :: f(nl,0:nx+1,0:ny+1,0:nz+1)


   ! define a tecplot object
   type(tecplot_time_file) :: plt_file
   integer,allocatable :: locations(:)
   integer,allocatable :: type_list(:)
   integer,allocatable :: shared_list(:)
   integer :: i,j,k,d
   real(kind=4), allocatable :: your_datas(:,:,:,:)
   real(kind=4) :: physics_time
   real :: xyz(3)


   real, allocatable :: blanking(:,:,:)
   allocate(blanking(nx,ny,nz))
   blanking=0
   where (lblanking) blanking=1


   allocate(your_datas(nx,ny,nz,num_of_variables))
   allocate(locations(num_of_variables))
   allocate(type_list(num_of_variables))
   allocate(shared_list(num_of_variables))

   ! locations = 0 means data in node, 1 means data in cell(not supported yet)
   locations = 0
   ! shared_list(i)=-1 means the i-th data is not shared in this zone. If shared_list(i)=m,
   ! it means the i-th data is shared with zone m in this file
   shared_list = -1
   ! type_list(i) = 0 means the i-th data is of type float. (Other data type not supported yet.)
   type_list = 1

   ! call init subroutine first
   ! nx, ny, nz means the dimension of the data
   ! 'x,y,z,u,v,w' is a string contains names of variables, must be divided by ','
   !call plt_file%init(filename,nx,ny,nz,'Tecplot File Title','x,y,z,u,v,w')
   call plt_file%init(filename,nx,ny,nz,'LBM3D output',trim(variables_string))

   ! for each zone, call the two subroutines
   ! physics_time can be any value, it will only be used when there are more than 1 zone in a file.
   call plt_file%write_zone_header(filename(4:9), physics_time, 0, locations)

   ! your_datas(:,:,:,1:3) =  x,y,z coordinates(Variable assignment is omitted in this example)
   ! ALL datas are stored in sequence like (((x(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
   ! set coordinate
   do d = 1, 3
      do concurrent(i=1:nx, j=1:ny, k=1:nz)
         xyz = [i-1., j-1., k-1.]
         your_datas(i,j,k,d) = xyz(d)
      end do
   end do

   ! set value
   your_datas(:,:,:,4)  = blanking(:,:,:)
   do d=1,27
      your_datas(:,:,:,4+d)  = f(:,:,:,d)
   enddo

   call plt_file%write_zone_data(type_list, shared_list, your_datas)

   ! before exit, you must call complete subroutine
   call plt_file%complete
   deallocate(blanking)

end subroutine

end module
