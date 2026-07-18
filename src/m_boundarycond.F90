! Revised boundary-condition driver and mixed-edge/corner completion.
!
! Assumptions:
!   ibnd/jbnd/kbnd = 0 periodic, 1 open inflow/outflow, 11:22 closed.
!   Open inflow/outflow is implemented only in i and j.
!   Closed face kernels operate on face interiors only.
!   Existing closed-closed edge kernels are retained.
!
! Required corrected existing modules:
!   m_boundary_i_inflow_kernel
!   m_boundary_j_inflow_kernel
!   m_boundary_i_periodic_kernel
!   m_boundary_j_periodic_kernel
!   m_boundary_k_periodic_kernel
!   m_boundary_i_closed_kernel
!   m_boundary_j_closed_kernel
!   m_boundary_k_closed_kernel
!   m_boundary_i_j_edges
!   m_boundary_i_k_edges
!   m_boundary_j_k_edges
!   m_boundary_closed_corners
!
module m_boundarycond
contains
subroutine boundarycond(f1,f2,uvel)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile, only : ibnd,jbnd,kbnd,rho0,udir
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime

   use m_boundary_i_inflow_kernel
   use m_boundary_j_inflow_kernel
   use m_boundary_ij_corner_kernel

   use m_boundary_i_periodic_kernel
   use m_boundary_j_periodic_kernel
   use m_boundary_k_periodic_kernel

   use m_boundary_i_closed_kernel
   use m_boundary_j_closed_kernel
   use m_boundary_k_closed_kernel

   use m_boundary_i_j_edges
   use m_boundary_i_k_edges
   use m_boundary_j_k_edges
   use m_boundary_closed_corners

   use m_boundary_mixed_edges
   use m_boundary_nongeneric_corners

   implicit none
   real, intent(inout) :: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: uvel(nz)

   real :: taperi(nx),taperj(ny),taperk(nz)
#ifdef _CUDA
   attributes(device) :: f1,f2,uvel,taperi,taperj,taperk
   integer :: tx,ty,tz,bx,by,bz
#endif
   integer, parameter :: icpu=11
   integer :: opt_i1,opt_iN,opt_j1,opt_jN,opt_k1,opt_kN
   integer :: opt_ij,opt_ik,opt_jk,opt_ijk

   call cpustart()

   if (kbnd == 1) error stop 'boundarycond: open boundary kbnd=1 is not implemented'

   call decode_closed(ibnd,opt_i1,opt_iN)
   call decode_closed(jbnd,opt_j1,opt_jN)
   call decode_closed(kbnd,opt_k1,opt_kN)

   taperi = 1.0
   taperj = 1.0
   taperk = 1.0

   !-----------------------------------------------------------------
   ! 1. Preliminary periodic sweep.
   !    This supplies valid tangential ghosts before closed-wall work.
   !-----------------------------------------------------------------
   call apply_periodic_sweep(f1)

   !-----------------------------------------------------------------
   ! 2. Open face interiors.
   !-----------------------------------------------------------------
   if (ibnd == 1) then
#ifdef _CUDA
      tx=1;   bx=1
      ty=8;   by=(ny+ty-1)/ty
      tz=8;   bz=(nz+tz-1)/tz
#endif
      call boundary_i_inflow_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,uvel,udir,taperj,taperk)
   endif

   if (jbnd == 1) then
#ifdef _CUDA
      tx=1;   bx=1
      ty=8;   by=(nx+ty-1)/ty
      tz=8;   bz=(nz+tz-1)/tz
#endif
      call boundary_j_inflow_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,uvel,udir,taperi,taperk)
   endif

   ! Open-open i-j corner lines, k=1:nz.
   if (ibnd == 1 .and. jbnd == 1) then
#ifdef _CUDA
      tx=32; bx=(nl+tx-1)/tx
      ty=4;  by=(nz+ty-1)/ty
#endif
      call boundary_ij_corner_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,1),dim3(tx,ty,1)>>>&
#endif
      &(f1)
   endif

   !-----------------------------------------------------------------
   ! 3. Closed face interiors.
   !-----------------------------------------------------------------
   if (opt_i1 /= 0) then
#ifdef _CUDA
      tx=1;   bx=1
      ty=nty; by=(ny+ty-1)/ty
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call boundary_i_closed_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,1,opt_i1)

      call boundary_i_closed_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,nx,opt_iN)
   endif

   if (opt_j1 /= 0) then
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=1;   by=1
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call boundary_j_closed_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,1,opt_j1)

      call boundary_j_closed_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,ny,opt_jN)
   endif

   if (opt_k1 /= 0) then
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=1;   bz=1
#endif
      call boundary_k_closed_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,1,opt_k1)

      call boundary_k_closed_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,nz,opt_kN)
   endif

   !-----------------------------------------------------------------
   ! 4. Closed-closed edges. The conservative rule is:
   !    free-slip only when both intersecting walls are free-slip;
   !    otherwise use no-slip at that edge.
   !-----------------------------------------------------------------
   opt_ij = merge(+1,-1,ibnd==22 .and. jbnd==22)
   opt_ik = merge(+1,-1,ibnd==22 .and. kbnd==22)
   opt_jk = merge(+1,-1,jbnd==22 .and. kbnd==22)

   if (ibnd>10 .and. jbnd>10) call boundary_i_j_edges(f1,f2,opt_ij)
   if (ibnd>10 .and. kbnd>10) call boundary_i_k_edges(f1,f2,opt_ik)
   if (jbnd>10 .and. kbnd>10) call boundary_j_k_edges(f1,f2,opt_jk)

   !-----------------------------------------------------------------
   ! 5. Mixed open-closed edges. Extrapolate the already completed
   !    closed edge along the open-boundary normal. This preserves the
   !    closed-wall reflection and avoids stale ghost populations.
   !-----------------------------------------------------------------
   call boundary_mixed_edges(f1,ibnd,jbnd,kbnd)

   !-----------------------------------------------------------------
   ! 6. Three-closed-wall corners.
   !-----------------------------------------------------------------
   if (ibnd>10 .and. jbnd>10 .and. kbnd>10) then
      opt_ijk = merge(+1,-1,ibnd==22 .and. jbnd==22 .and. kbnd==22)
      call boundary_closed_corners(f1,f2,opt_ijk)
   else
      call boundary_nongeneric_corners(f1,ibnd,jbnd,kbnd)
   endif

   !-----------------------------------------------------------------
   ! 7. Final periodic sweeps. Two passes are intentional: one pass
   !    can update an edge from a still-old corner in another periodic
   !    direction; the second closes all periodic intersections.
   !-----------------------------------------------------------------
   call apply_periodic_sweep(f1)
   call apply_periodic_sweep(f1)

   call cpufinish(icpu)

contains

   subroutine decode_closed(bnd,opt1,optN)
      integer, intent(in)  :: bnd
      integer, intent(out) :: opt1,optN
      select case(bnd)
      case(11); opt1=-1; optN=-1
      case(12); opt1=-1; optN=+1
      case(21); opt1=+1; optN=-1
      case(22); opt1=+1; optN=+1
      case default; opt1=0; optN=0
      end select
   end subroutine decode_closed

   subroutine apply_periodic_sweep(f)
      real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
      attributes(device) :: f
#endif
      if (ibnd==0) then
#ifdef _CUDA
         tx=1;  bx=1
         ty=32; by=(ny+2+ty-1)/ty
         tz=8;  bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
         &(f,nl)
      endif
      if (jbnd==0) then
#ifdef _CUDA
         tx=256; bx=(nx+2+tx-1)/tx
         ty=1;   by=1
         tz=1;   bz=(nz+2+tz-1)/tz
#endif
         call boundary_j_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
         &(f,nl)
      endif
      if (kbnd==0) then
#ifdef _CUDA
         tx=256; bx=(nl*(nx+2)+tx-1)/tx
         ty=1;   by=(ny+2+ty-1)/ty
         tz=1;   bz=1
#endif
         call boundary_k_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
         &(f,nl)
      endif
   end subroutine apply_periodic_sweep

end subroutine boundarycond
end module m_boundarycond


