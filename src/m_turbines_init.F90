module m_turbines_init
   integer, save :: iradius
   integer, parameter, public :: ieps=5            ! number of gridcells for smoothing actuatorline forcing in i-dir
   real, save :: theta=0.0      ! start angle for turbine rotation

   real, allocatable  :: turbine_df(:,:,:,:,:)      ! Turbine forcing used in applyturbines
   real, allocatable  :: dfeq1(:,:,:,:)
   real, allocatable  :: dfeq2(:,:,:,:)
   real, allocatable  :: force(:,:,:,:)     ! work array for computing the turbine force
   real, allocatable  :: du(:,:,:)          ! turbine forced u velocity
   real, allocatable  :: dv(:,:,:)          ! turbine forced v velocity
   real, allocatable  :: dw(:,:,:)          ! turbine forced w velocity
   real, allocatable  :: vel(:,:,:,:)
   real, allocatable  :: rtmp(:,:,:)
   real, allocatable  :: cx(:), cy(:), cz(:)
   real, allocatable  :: slice_u(:,:)
   real, allocatable  :: slice_v(:,:)
   real, allocatable  :: slice_w(:,:)
   real, allocatable  :: slice_r(:,:)
#ifdef _CUDA
   attributes(device) :: turbine_df
   attributes(device) :: dfeq1
   attributes(device) :: dfeq2
   attributes(device) :: force
   attributes(device) :: du
   attributes(device) :: dv
   attributes(device) :: dw
   attributes(device) :: vel
   attributes(device) :: rtmp
   attributes(device) :: cx, cy, cz
   attributes(device) :: slice_u
   attributes(device) :: slice_v
   attributes(device) :: slice_w
   attributes(device) :: slice_r
#endif

contains
subroutine turbines_init
   use mod_dimensions
   use m_readinfile, only : nturbines
   use mod_D3Q27setup
   implicit none
   integer l
   if (.not. allocated(turbine_df))  allocate(turbine_df(1:nl,-ieps:ieps,1:ny,1:nz,nturbines))
   if (.not. allocated(dfeq1))       allocate(dfeq1(nl,-ieps:ieps,ny,nz)            )
   if (.not. allocated(dfeq2))       allocate(dfeq2(nl,-ieps:ieps,ny,nz)            )
   if (.not. allocated(force))       allocate(force(0:ieps,ny,nz,3))
   if (.not. allocated(du))          allocate(du(-ieps:ieps,ny,nz) )
   if (.not. allocated(dv))          allocate(dv(-ieps:ieps,ny,nz) )
   if (.not. allocated(dw))          allocate(dw(-ieps:ieps,ny,nz) )
   if (.not. allocated(vel))         allocate(vel(3,-ieps:ieps,ny,nz) )
   if (.not. allocated(rtmp))        allocate(rtmp(-ieps:ieps,ny,nz) )
   if (.not. allocated(cx))          allocate(cx(nl))
   if (.not. allocated(cy))          allocate(cy(nl))
   if (.not. allocated(cz))          allocate(cz(nl))
   if (.not. allocated(slice_u))     allocate(slice_u(ny,nz))
   if (.not. allocated(slice_v))     allocate(slice_v(ny,nz))
   if (.not. allocated(slice_w))     allocate(slice_w(ny,nz))
   if (.not. allocated(slice_r))     allocate(slice_r(ny,nz))

   do l=1,nl
     cx(l)=real(cxs_h(l))
     cy(l)=real(cys_h(l))
     cz(l)=real(czs_h(l))
   enddo
end subroutine
end module


