module m_tecout
implicit none

contains
subroutine tecout(filetype,filename,it,variables_string,num_of_variables,lblanking,rho,u,v,w,Ti)
   use mod_dimensions
   use m_tecplot
   use m_readinfile
#ifdef MPI
   use m_mpi_decomp_init, only : j_start
#endif
   implicit none
   character(len=*), intent(in) :: filename         ! Output filename
   integer,          intent(in) :: filetype           ! Save only geometry when filetype=1
   integer,          intent(in) :: it               ! Timestep index
   character(len=*), intent(in) :: variables_string ! Variables to print separated by ,
   integer,          intent(in) :: num_of_variables ! Number of vaiables to print

   real,     intent(in)             :: rho(0:nx+1,0:ny+1,0:nz+1)       ! fluid density
   real,     intent(in)             :: u(0:nx+1,0:ny+1,0:nz+1)         ! x component of fluid velocity
   real,     intent(in)             :: v(0:nx+1,0:ny+1,0:nz+1)         ! y component of fluid velocity
   real,     intent(in)             :: w(0:nx+1,0:ny+1,0:nz+1)         ! z component of fluid velocity
   logical,  intent(in)             :: lblanking(0:nx+1,0:ny+1,0:nz+1) ! blanking
   real,     intent(in), optional   :: Ti(0:nx+1,0:ny+1,0:nz+1)        ! Turbulent kinetic enery


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
#ifndef MPI
   integer :: j_start=1
#endif

   physics_time=real(it)
   print '(5a,f10.2)','tecout: ',trim(filename),' ',trim(variables_string),' iteration=',physics_time

   if ((filetype == 0) .or. (filetype == 1)) then
      allocate(blanking(nx,ny+1,nz))
      blanking=0.0
      do k=1,nz
      do j=1,ny+1
      do i=1,nx
         if (lblanking(i,j,k)) blanking(i,j,k)=1.0
      enddo
      enddo
      enddo
   endif

   allocate(your_datas(nx,ny+1,nz,num_of_variables))
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
   call plt_file%init(filename,nx,ny+1,nz,filetype,'LBM3D output',filetype,trim(variables_string))

   ! for each zone, call the two subroutines
   ! physics_time can be any value, it will only be used when there are more than 1 zone in a file.
   call plt_file%write_zone_header(filename(4:9), physics_time, 0, locations)

   ! your_datas(:,:,:,1:3) =  x,y,z coordinates(Variable assignment is omitted in this example)
   ! ALL datas are stored in sequence like (((x(ix,iy,iz),ix=1,nx),iy=1,ny),iz=1,nz)
   ! set coordinate

! static grid data
   if ((filetype == 0) .or. (filetype == 1)) then
      do k = 1, nz
         do j = 1, ny+1
            do i = 1, nx
               your_datas(i,j,k,1) = real(i - 1)
               your_datas(i,j,k,2) = real(j_start + j - 1)
               your_datas(i,j,k,3) = real(k - 1)
             ! your_datas(i,j,k,1) = (i - 1) * dx
             ! your_datas(i,j,k,2) = (jstart + j - 1) * dy
             ! your_datas(i,j,k,3) = (k - 1) * dz
            end do
         end do
      end do

      your_datas(:,:,:,4) = blanking(:,:,:)
   endif

! solution data
   if ((filetype == 0) .or. (filetype == 2)) then
      if (filetype == 2) dd=0
      dd=dd+1; your_datas(1:nx,1:ny+1,1:nz,dd)  = rho(1:nx,1:ny+1,1:nz)
      dd=dd+1; your_datas(1:nx,1:ny+1,1:nz,dd)  =   u(1:nx,1:ny+1,1:nz)
      dd=dd+1; your_datas(1:nx,1:ny+1,1:nz,dd)  =   v(1:nx,1:ny+1,1:nz)
      dd=dd+1; your_datas(1:nx,1:ny+1,1:nz,dd)  =   w(1:nx,1:ny+1,1:nz)
      if (present(Ti)) then
         dd=dd+1; your_datas(1:nx,1:ny+1,1:nz,dd)  = Ti(1:nx,1:ny+1,1:nz)
      endif
   endif

   call plt_file%write_zone_data(type_list, shared_list, your_datas)

   ! before exit, you must call complete subroutine
   call plt_file%complete

   if (allocated(blanking)) deallocate(blanking)

end subroutine

end module
