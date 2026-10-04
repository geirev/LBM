module m_inflow_turbulence_init

   integer, parameter :: iturb_pos=1
   integer :: nturb_perim=0

   integer, allocatable :: iturb_i(:), iturb_j(:), iturb_s(:), iturb_face(:)
   real,    allocatable :: uu(:,:,:), vv(:,:,:), ww(:,:,:), rr(:,:,:)

#ifdef _CUDA
   attributes(device) :: iturb_i,iturb_j,iturb_s,iturb_face
   attributes(device) :: uu,vv,ww,rr
#endif

contains

subroutine inflow_turbulence_init

   use mod_dimensions, only : nx,ny,nz,nyg
   use m_readinfile,   only : nrturb
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif

   implicit none

   integer :: i0,i1,j0g,j1g
   integer :: i,jg,jloc,s
   integer :: nperim,np
   integer, allocatable :: i_h(:),j_h(:),s_h(:),face_h(:)
   integer :: joff

   if (nrturb <= 0) stop 'Call read_infile before inflow_turbulence_init'

   ! Forcing rectangle, iturb_pos cells inside the global open boundaries.
   i0  = 1+iturb_pos
   i1  = nx-iturb_pos
   j0g = 1+iturb_pos
   j1g = nyg-iturb_pos

   if (i1 <= i0 .or. j1g <= j0g) &
      stop 'inflow_turbulence_init: forcing rectangle is too small'

   ! Number of unique points around the global perimeter.
   nperim = 2*(i1-i0+1) + 2*(j1g-j0g+1) - 4

#ifdef MPI
   joff = mpi_rank*ny
#else
   joff = 0
#endif

   ! First pass: count perimeter points owned by this j tile.
   nturb_perim=0
   s=0

   ! Bottom: left -> right, including both corners.
   jg=j0g
   do i=i0,i1
      s=s+1
      if (jg > joff .and. jg <= joff+ny) nturb_perim=nturb_perim+1
   enddo

   ! Right: bottom -> top, excluding bottom corner.
   i=i1
   do jg=j0g+1,j1g
      s=s+1
      if (jg > joff .and. jg <= joff+ny) nturb_perim=nturb_perim+1
   enddo

   ! Top: right -> left, excluding right corner.
   jg=j1g
   do i=i1-1,i0,-1
      s=s+1
      if (jg > joff .and. jg <= joff+ny) nturb_perim=nturb_perim+1
   enddo

   ! Left: top -> bottom, excluding both corners.
   i=i0
   do jg=j1g-1,j0g+1,-1
      s=s+1
      if (jg > joff .and. jg <= joff+ny) nturb_perim=nturb_perim+1
   enddo

   if (s /= nperim) stop 'inflow_turbulence_init: perimeter count error'

   allocate(iturb_i(nturb_perim),iturb_j(nturb_perim))
   allocate(iturb_s(nturb_perim),iturb_face(nturb_perim))
   allocate(i_h(nturb_perim),j_h(nturb_perim),s_h(nturb_perim),face_h(nturb_perim))
   allocate(uu(nturb_perim,nz,0:nrturb))
   allocate(vv(nturb_perim,nz,0:nrturb))
   allocate(ww(nturb_perim,nz,0:nrturb))
   allocate(rr(nturb_perim,nz,0:nrturb))

   ! Second pass: local mapping.  Face identifiers are
   ! 1=bottom(y-min), 2=right(x-max), 3=top(y-max), 4=left(x-min).
   nturb_perim=0
   s=0

   jg=j0g
   do i=i0,i1
      s=s+1
      if (jg > joff .and. jg <= joff+ny) then
         nturb_perim=nturb_perim+1
         jloc=jg-joff
         i_h(nturb_perim)=i
         j_h(nturb_perim)=jloc
         s_h(nturb_perim)=s
         face_h(nturb_perim)=1
      endif
   enddo

   i=i1
   do jg=j0g+1,j1g
      s=s+1
      if (jg > joff .and. jg <= joff+ny) then
         nturb_perim=nturb_perim+1
         jloc=jg-joff
         i_h(nturb_perim)=i
         j_h(nturb_perim)=jloc
         s_h(nturb_perim)=s
         face_h(nturb_perim)=2
      endif
   enddo

   jg=j1g
   do i=i1-1,i0,-1
      s=s+1
      if (jg > joff .and. jg <= joff+ny) then
         nturb_perim=nturb_perim+1
         jloc=jg-joff
         i_h(nturb_perim)=i
         j_h(nturb_perim)=jloc
         s_h(nturb_perim)=s
         face_h(nturb_perim)=3
      endif
   enddo

   i=i0
   do jg=j1g-1,j0g+1,-1
      s=s+1
      if (jg > joff .and. jg <= joff+ny) then
         nturb_perim=nturb_perim+1
         jloc=jg-joff
         i_h(nturb_perim)=i
         j_h(nturb_perim)=jloc
         s_h(nturb_perim)=s
         face_h(nturb_perim)=4
      endif
   enddo

   iturb_i=i_h
   iturb_j=j_h
   iturb_s=s_h
   iturb_face=face_h
   deallocate(i_h,j_h,s_h,face_h)

   uu=0.0; vv=0.0; ww=0.0; rr=0.0

end subroutine inflow_turbulence_init

end module m_inflow_turbulence_init
