module m_tecout
implicit none

contains
subroutine tecout(filetype,filename,it,variables_string,num_of_variables,lblanking,rho,u,v,w,Ti)
   use mod_dimensions
   use m_tecplot
   use m_readinfile
   implicit none
   character(len=*), intent(in) :: filename         ! Output filename
   integer,          intent(in) :: filetype           ! Save only geometry when filetype=1
   integer,          intent(in) :: it               ! Timestep index
   character(len=*), intent(in) :: variables_string ! Variables to print separated by ,
   integer,          intent(in) :: num_of_variables ! Number of vaiables to print

   real,     intent(in)             :: rho(nx,ny,nz)       ! fluid density
   real,     intent(in)             :: u(nx,ny,nz)         ! x component of fluid velocity
   real,     intent(in)             :: v(nx,ny,nz)         ! y component of fluid velocity
   real,     intent(in)             :: w(nx,ny,nz)         ! z component of fluid velocity
   logical,  intent(in)             :: lblanking(0:nx+1,0:ny+1,0:nz+1) ! blanking
   real,     intent(in), optional   :: Ti(nx,ny,nz)        ! Turbulent kinetic enery


   ! define a tecplot object
   type(tecplot_time_file) :: plt_file
   integer,allocatable :: locations(:)
   integer,allocatable :: type_list(:)
   integer,allocatable :: shared_list(:)
   integer :: i,j,k,d
   real(kind=4), allocatable :: your_datas(:,:,:,:)
   real(kind=4) :: physics_time
   real :: xyz(3)
   integer dd
   real, allocatable :: blanking(:,:,:)

   physics_time=real(it)
   print '(5a,f10.2)','tecout: ',trim(filename),' ',trim(variables_string),' iteration=',physics_time

   if ((filetype == 0) .or. (filetype == 1)) then
      allocate(blanking(nx,ny,nz))
      blanking=0.0
      do k=1,nz
      do j=1,ny
      do i=1,nx
         if (lblanking(i,j,k)) blanking(i,j,k)=1.0
      enddo
      enddo
      enddo
   endif

   allocate(your_datas(nx,ny,nz,num_of_variables))
   allocate(locations(num_of_variables))
   allocate(type_list(num_of_variables))
   allocate(shared_list(num_of_variables))

   ! locations = 0 means data in node, 1 means data in cell(not supported yet)
   locations = 0
   ! shared_list(i)=-1 means the i-th data is not shared in this zone. If shared_list(i)=m,
   ! it means the i-th data is shared with zone m in this file
      shared_list = -1
   ! type_list(i) = 1 means the i-th data is of type float. (Other data type not supported yet.)
   type_list = 1

   ! call init subroutine first
   ! nx, ny, nz means the dimension of the data
   ! 'x,y,z,u,v,w' is a string contains names of variables, must be divided by ','
   call plt_file%init(filename,nx,ny,nz,filetype,'LBM3D output',filetype,trim(variables_string))

   ! for each zone, call the two subroutines
   ! physics_time can be any value, it will only be used when there are more than 1 zone in a file.
   call plt_file%write_zone_header(filename(4:9), physics_time, 0, locations)

   ! your_datas(:,:,:,1:3) =  x,y,z coordinates(Variable assignment is omitted in this example)
   ! ALL datas are stored in sequence like (((x(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
   ! set coordinate

! static grid data
   if ((filetype == 0) .or. (filetype == 1)) then
      do d = 1, 3
         do concurrent(i=1:nx, j=1:ny, k=1:nz)
            xyz = [i-1., j-1., k-1.]
            your_datas(i,j,k,d) = xyz(d)
         end do
      end do
      dd=3+1; your_datas(:,:,:,dd)  = blanking(:,:,:)
   endif

! solution data
   if ((filetype == 0) .or. (filetype == 2)) then
      if (filetype == 2) dd=0
      dd=dd+1; your_datas(:,:,:,dd)  = rho(:,:,:)
      dd=dd+1; your_datas(:,:,:,dd)  = u(:,:,:)
      dd=dd+1; your_datas(:,:,:,dd)  = v(:,:,:)
      dd=dd+1; your_datas(:,:,:,dd)  = w(:,:,:)
      if (present(Ti)) then
         dd=dd+1; your_datas(:,:,:,dd)  = Ti(:,:,:)
      endif
   endif

   call plt_file%write_zone_data(type_list, shared_list, your_datas)

   ! before exit, you must call complete subroutine
   call plt_file%complete

   if (allocated(blanking)) deallocate(blanking)

end subroutine

end module
