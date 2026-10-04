module m_inflow_turbulence_init

   integer, parameter :: iturb_pos=10
   integer, parameter :: iturb_radius=2

   real, allocatable  :: cx(:), cy(:), cz(:)
   real, allocatable  :: uu(:,:,:), vv(:,:,:), ww(:,:,:), rr(:,:,:)

#ifdef _CUDA
   attributes(device) :: cx, cy, cz, uu, vv, ww, rr
#endif

contains
subroutine inflow_turbulence_init
   use mod_dimensions, only : ny,nz
   use mod_D3Q27setup, only : nl, cxs_h, cys_h, czs_h
   use m_readinfile,   only : nrturb
   implicit none
   integer :: l

   if (nrturb <= 0) stop 'Call read_infile before inflow_turbulence_init'

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

