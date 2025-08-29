module m_inflow_turbulence_init
   integer, parameter :: iturb_pos=10
   integer, parameter :: iturb_radius=2

   real, allocatable  :: turbulence_df(:,:,:)

! Persistent local device arrays
   real, allocatable  :: vel(:,:,:,:), rtmp(:,:,:)
   real, allocatable  :: dfeq1(:,:,:,:), dfeq2(:,:,:,:)
   real, allocatable  :: cx(:), cy(:), cz(:)

#ifdef _CUDA
   attributes(device) :: turbulence_df
   attributes(device) :: vel, rtmp
   attributes(device) :: dfeq1, dfeq2
   attributes(device) :: cx, cy, cz
#endif

! Stochastic input field on inflow boundary
   real, allocatable :: uu(:,:,:)
   real, allocatable :: vv(:,:,:)
   real, allocatable :: ww(:,:,:)
   real, allocatable :: rr(:,:,:)
#ifdef _CUDA
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
   attributes(device) :: rr
#endif


contains

subroutine inflow_turbulence_init
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile, only : nrturb
   implicit none
   integer l
   if (nrturb==0) then
      print *,'make sure read_infile is called before inflow_turbulence_init and nrturb>0'
      stop
   endif

   if (.not. allocated(turbulence_df)) allocate(turbulence_df(1:nl,1:ny,1:nz))
   if (.not. allocated(vel))           allocate(vel(3,1,ny,nz))
   if (.not. allocated(rtmp))          allocate(rtmp(1,ny,nz))
   if (.not. allocated(dfeq1))         allocate(dfeq1(nl,1,ny,nz))
   if (.not. allocated(dfeq2))         allocate(dfeq2(nl,1,ny,nz))
   if (.not. allocated(cx))            allocate(cx(nl))
   if (.not. allocated(cy))            allocate(cy(nl))
   if (.not. allocated(cz))            allocate(cz(nl))

   if (.not. allocated(uu))            allocate(uu(ny,nz,0:nrturb))
   if (.not. allocated(vv))            allocate(vv(ny,nz,0:nrturb))
   if (.not. allocated(ww))            allocate(ww(ny,nz,0:nrturb))
   if (.not. allocated(rr))            allocate(rr(ny,nz,0:nrturb))

   do l=1,nl
     cx(l)=real(cxs_h(l))
     cy(l)=real(cys_h(l))
     cz(l)=real(czs_h(l))
   enddo
end subroutine
end module


