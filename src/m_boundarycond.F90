module m_boundarycond
! Important: periodic BCs in j-direction must be applied
! before AND after closed BCs in k-direction.
!
! First application ensures that j-ghost cells are correct
! before bounce-back in k-direction.
!
! Second application ensures that the updated k=0,k=nz+1
! boundary populations are correctly wrapped in j-direction,
!
! maintaining consistency at edges and corners.
! Same logic applies if periodic in ks-direction  and closed in j-direction.
contains
subroutine boundarycond(f1,f2,uvel,tracer,pottemp)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile,   only : ibnd,jbnd,kbnd,rho0,udir,iablvisc,istable
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime

   use m_boundary_i_inflow_kernel
   use m_boundary_i_inflow_edges

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

   implicit none
   real, intent(inout):: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout):: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout):: tracer(ntracer,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout):: pottemp(0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: uvel(nz)

   real, taperj(ny),taperk(nz)!, real, dist, x, width
   real, parameter   :: pi=3.1415927410125732
#ifdef _CUDA
   attributes(device) :: f1
   attributes(device) :: f2
   attributes(device) :: tracer
   attributes(device) :: pottemp
   attributes(device) :: uvel
   attributes(device) :: taperj,taperk
#endif
   integer, parameter :: icpu=11
   integer :: opt_i1, opt_j1, opt_k1
   integer :: opt_iN, opt_jN, opt_kN
   !integer j,k
#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif
   
   call cpustart()


! Update edges for inflow conditions for periodic boundary conditions in i-direction
   if (ibnd==0) then
#ifdef _CUDA
      tx=1 ; bx=1
      ty=32; by=(ny+2+ty-1)/ty
      tz=8 ; bz=(nz+2+tz-1)/tz
#endif
      call boundary_i_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1,nl)

      if (ntracer > 0) then
#ifdef _CUDA
         tx=1;  bx=1
         ty=32; by=(ny+2+ty-1)/ty
         tz=8;  bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(tracer,ntracer)
      endif

      if (iablvisc > 0) then
#ifdef _CUDA
         tx=1  ; bx=1
         ty=32 ; by=(ny+2+ty-1)/ty
         tz=8  ; bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
      endif
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
#ifdef _CUDA
      tx=256; bx=(nx+2+tx-1)/tx
      ty=1;   by=1
      tz=1;   bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1,nl)

      if (ntracer > 0) then
#ifdef _CUDA
      tx=256; bx=(nx+2+tx-1)/tx
      ty=1;   by=1
      tz=1;   bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(tracer,ntracer)
      endif

      if (iablvisc > 0) then
#ifdef _CUDA
         tx=256; bx=(nx+2+tx-1)/tx
         ty=1;   by=1
         tz=1;   bz=(nz+2+tz-1)/tz
#endif
         call boundary_j_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
      endif
   endif

   if (kbnd==0) then
#ifdef _CUDA
      tx=256; bx=(nl*(nx+2)+tx-1)/tx
      ty=1;   by=(ny+2+ty-1)/ty
      tz=1;   bz=1
#endif
      call boundary_k_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1)
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Setting up for closed boundaries
! i-direction
   select case(ibnd)
   case(11)     ; opt_i1 = -1; opt_iN = -1
   case(12)     ; opt_i1 = -1; opt_iN = +1
   case(21)     ; opt_i1 = +1; opt_iN = -1
   case(22)     ; opt_i1 = +1; opt_iN = +1
   case default ; opt_i1 =  0; opt_iN =  0
   end select

   if (opt_i1 /= 0) then
#ifdef _CUDA
      tx=1; bx=1
      ty=nty;  by=(ny+ty-1)/ty
      tz=ntz;  bz=(nz+tz-1)/tz
#endif

      call boundary_i_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,1,opt_i1)

      call boundary_i_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,nx,opt_iN)
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! j-direction
   select case(jbnd)
   case(11)     ; opt_j1 = -1; opt_jN = -1
   case(12)     ; opt_j1 = -1; opt_jN = +1
   case(21)     ; opt_j1 = +1; opt_jN = -1
   case(22)     ; opt_j1 = +1; opt_jN = +1
   case default ; opt_j1 =  0; opt_jN =  0
   end select

   if (opt_j1 /= 0) then
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+tz-1)/tz
#endif

      call boundary_j_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,1,opt_j1)

      call boundary_j_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,ny,opt_jN)
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! k-direction
   select case(kbnd)
   case(11)     ; opt_k1 = -1; opt_kN = -1
   case(12)     ; opt_k1 = -1; opt_kN = +1
   case(21)     ; opt_k1 = +1; opt_kN = -1
   case(22)     ; opt_k1 = +1; opt_kN = +1
   case default ; opt_k1 =  0; opt_kN =  0
   end select

   if (opt_k1 /= 0) then
! No-slip bounce-back condition for k=1
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=1; bz=1
#endif

      call boundary_k_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,1,opt_k1)

      call boundary_k_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,nz,opt_kN)
   endif


!======================================================================
!  Edges and corners for closed boundaries (i,j,k)
!
!  Rule:
!     - Only process edges if both directions are closed ( > 10 )
!     - Free-slip only if *both* directions are free-slip (== 22)
!     - Otherwise no-slip
!======================================================================

!-------------------------------------------------------------
! 1. j–k edges and jk-corners
!-------------------------------------------------------------
   if (jbnd>10 .and. kbnd>10) then

      if (jbnd==22 .and. kbnd==22) then
         print *, "JK edges: free-slip"
         call boundary_j_k_edges(f1,f2,1)
      else
         print *, "JK edges: no-slip"
         call boundary_j_k_edges(f1,f2,-1)
      endif

   endif


!-------------------------------------------------------------
! 2. i–j edges and ij-corners
!-------------------------------------------------------------
   if (ibnd>10 .and. jbnd>10) then

      if (ibnd==22 .and. jbnd==22) then
         print *, "IJ edges: free-slip"
         call boundary_i_j_edges(f1,f2,1)        ! NEW routine
      else
         print *, "IJ edges: no-slip"
         call boundary_i_j_edges(f1,f2,-1)       ! NEW routine
      endif

   endif


!-------------------------------------------------------------
! 3. i–k edges and ik-corners
!-------------------------------------------------------------
   if (ibnd>10 .and. kbnd>10) then

      if (ibnd==22 .and. kbnd==22) then
         print *, "IK edges: free-slip"
         call boundary_i_k_edges(f1,f2,1)        ! NEW routine
      else
         print *, "IK edges: no-slip"
         call boundary_i_k_edges(f1,f2,-1)       ! NEW routine
      endif

   endif

   if (ibnd>10 .and. jbnd>10 .and. kbnd > 10) then
      print *, "closed corners"
      call boundary_closed_corners(f1,f2,-1)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! i-inflow boundary condition
   if (ibnd==1) then
!      width=4.0
!      do j=1,ny
!         dist=min(abs(j-1),abs(ny-j))
!         if (dist >= width) then
!            taperj(j) = 1.0
!         else
!            x = dist / width
!            taperj(j) = cos(0.5*pi*x)**2
!         endif
!      enddo
!      do k=1,nz
!         dist=min(abs(k-1),abs(nz-k))
!         if (dist >= width) then
!            taperk(k) = 1.0
!         else
!            x = dist / width
!            taperk(k) = cos(0.5*pi*x)**2
!         endif
!      enddo
      taperj=1.0
      taperk=1.0


#ifdef _CUDA
      tx=1;   bx=1
      ty=8;   by=(ny+ty-1)/ty
      tz=8;   bz=(nz+tz-1)/tz
#endif
      call boundary_i_inflow_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1,uvel,rho0,udir,tracer,pottemp,iablvisc,jbnd,kbnd,taperj,taperk)

! Periodic bnd cond for pottemp in unstable case
      if (iablvisc > 0 .and. istable == -1) then
#ifdef _CUDA
         tx=1  ; bx=1
         ty=32 ; by=(ny+2+ty-1)/ty
         tz=8  ; bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
      endif

! Update edges for inflow conditions
    !  call boundary_i_inflow_edges(f1,f2)
   endif






! Update edges for inflow conditions for periodic boundary conditions in i-direction
   if (ibnd==0) then
#ifdef _CUDA
      tx=1 ; bx=1
      ty=32; by=(ny+2+ty-1)/ty
      tz=8 ; bz=(nz+2+tz-1)/tz
#endif
      call boundary_i_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1,nl)

      if (ntracer > 0) then
#ifdef _CUDA
         tx=1;  bx=1
         ty=32; by=(ny+2+ty-1)/ty
         tz=8;  bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(tracer,ntracer)
      endif

      if (iablvisc > 0) then
#ifdef _CUDA
         tx=1  ; bx=1
         ty=32 ; by=(ny+2+ty-1)/ty
         tz=8  ; bz=(nz+2+tz-1)/tz
#endif
         call boundary_i_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
      endif
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
#ifdef _CUDA
      tx=256; bx=(nx+2+tx-1)/tx
      ty=1;   by=1
      tz=1;   bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1,nl)

      if (ntracer > 0) then
#ifdef _CUDA
      tx=256; bx=(nx+2+tx-1)/tx
      ty=1;   by=1
      tz=1;   bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(tracer,ntracer)
      endif

      if (iablvisc > 0) then
#ifdef _CUDA
         tx=256; bx=(nx+2+tx-1)/tx
         ty=1;   by=1
         tz=1;   bz=(nz+2+tz-1)/tz
#endif
         call boundary_j_periodic_kernel&
#ifdef _CUDA
         &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
         &(pottemp,1)
      endif
   endif

   if (kbnd==0) then
#ifdef _CUDA
      tx=256; bx=(nl*(nx+2)+tx-1)/tx
      ty=1;   by=(ny+2+ty-1)/ty
      tz=1;   bz=1
#endif
      call boundary_k_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1)
   endif


   call cpufinish(icpu)

end subroutine
end module
