! Revised boundary-condition driver and mixed-edge/corner completion.
!
! Assumptions:
!   ibnd/jbnd/kbnd = 0 periodic, 1 open inflow/outflow, 11:22 closed.
!   Open inflow/outflow is implemented only in i and j.
!   Closed face kernels operate on face interiors only.
!
!
module m_boundarycond
contains
subroutine boundarycond(f1,f2,rho,uvel)
   use mod_dimensions
   use mod_D3Q27setup, only : nl,cs2,cs4,cs6
   use m_readinfile, only : ibnd,jbnd,kbnd,udir,rho0,ibgk
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
#ifdef MPI
   use mpi
   use m_mpi_decomp_init, only : north, south
#endif
   use m_wtime

   use m_boundary_inflow_i_kernel
   use m_boundary_inflow_j_kernel
   use m_boundary_inflow_edges_ij_kernel

   use m_boundary_periodic_i_kernel
   use m_boundary_periodic_j_kernel
   use m_boundary_periodic_k_kernel

   use m_boundary_closed_i_kernel
   use m_boundary_closed_j_kernel
   use m_boundary_closed_k_kernel

   use m_boundary_edges_i_j
   use m_boundary_edges_i_k
   use m_boundary_edges_j_k
   use m_boundary_edges_mixed

   use m_boundary_corner_closed
   use m_boundary_corner_nongeneric

   implicit none
   real, intent(inout) :: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout) :: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: rho(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: uvel(nz)

#ifdef _CUDA
   attributes(device) :: f1,f2,uvel,rho
   integer :: tx,ty,tz,bx,by,bz
#endif
   integer, parameter :: icpu=11
   integer :: opt_i1,opt_iN,opt_j1,opt_jN,opt_k1,opt_kN
   integer :: opt_ij,opt_ik,opt_jk,opt_ijk
   real, parameter :: pi=acos(-1.0)
   real, parameter :: rho_relax=0.5
   real, parameter :: inv1cs2 = 1.0/(cs2)
   real, parameter :: inv2cs4 = 1.0/(2.0*cs4)
   real, parameter :: inv6cs6 = 1.0/(6.0*cs6)
   logical j0_is_phys, jN_is_phys

   call cpustart()

   if (kbnd == 1) error stop 'boundarycond: open boundary kbnd=1 is not implemented'

   call decode_closed(ibnd,opt_i1,opt_iN)
   call decode_closed(jbnd,opt_j1,opt_jN)
   call decode_closed(kbnd,opt_k1,opt_kN)

#ifdef MPI
      j0_is_phys = (south == MPI_PROC_NULL)
      jN_is_phys = (north == MPI_PROC_NULL)
#else
      j0_is_phys = .true.
      jN_is_phys = .true.
#endif

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
      call boundary_inflow_i_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,uvel,udir,rho0,rho_relax,inv1cs2,inv2cs4,inv6cs6,ibgk)
   endif


   if (jbnd == 1) then
#ifdef _CUDA
      tx=1;   bx=1
      ty=8;   by=(nx+ty-1)/ty
      tz=8;   bz=(nz+tz-1)/tz
#endif
      call boundary_inflow_j_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,uvel,udir,rho0,rho_relax,inv1cs2,inv2cs4,inv6cs6,ibgk, &
        j0_is_phys,jN_is_phys)
   endif

   ! Open-open i-j corner lines, k=1:nz.
   if (ibnd == 1 .and. jbnd == 1) then
#ifdef _CUDA
      tx=32; bx=(nl+tx-1)/tx
      ty=4;  by=(nz+ty-1)/ty
#endif
      call boundary_inflow_edges_ij_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,1),dim3(tx,ty,1)>>>&
#endif
      &(f1,uvel,udir,rho0,rho_relax,inv1cs2,inv6cs6,ibgk,j0_is_phys,jN_is_phys)
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
      call boundary_closed_i_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,1,opt_i1)

      call boundary_closed_i_kernel&
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
      call boundary_closed_j_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,1,opt_j1)

      call boundary_closed_j_kernel&
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
      call boundary_closed_k_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
      &(f1,f2,1,opt_k1)

      call boundary_closed_k_kernel&
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

   if (ibnd>10 .and. jbnd>10) call boundary_edges_i_j(f1,f2,opt_ij)
   if (ibnd>10 .and. kbnd>10) call boundary_edges_i_k(f1,f2,opt_ik)
   if (jbnd>10 .and. kbnd>10) call boundary_edges_j_k(f1,f2,opt_jk)

   !-----------------------------------------------------------------
   ! 5. Mixed open-closed edges. Extrapolate the already completed
   !    closed edge along the open-boundary normal. This preserves the
   !    closed-wall reflection and avoids stale ghost populations.
   !-----------------------------------------------------------------
   call boundary_edges_mixed(f1,ibnd,jbnd,kbnd,j0_is_phys,jN_is_phys)

   !-----------------------------------------------------------------
   ! 6. Three-closed-wall corners.
   !-----------------------------------------------------------------
   if (ibnd>10 .and. jbnd>10 .and. kbnd>10) then
      opt_ijk = merge(+1,-1,ibnd==22 .and. jbnd==22 .and. kbnd==22)
      call boundary_corner_closed(f1,f2,opt_ijk)
   else
      call boundary_corner_nongeneric(f1,uvel,udir,rho0,rho_relax,ibnd,jbnd,kbnd,j0_is_phys,jN_is_phys)
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
         call boundary_periodic_i_kernel&
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
         call boundary_periodic_j_kernel&
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
         call boundary_periodic_k_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz),dim3(tx,ty,tz)>>>&
#endif
         &(f,nl)
      endif
   end subroutine apply_periodic_sweep

end subroutine boundarycond
end module m_boundarycond


