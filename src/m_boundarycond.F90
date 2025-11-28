module m_boundarycond
contains
subroutine boundarycond(f1,f2,uvel,tracer)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile,   only : ibnd,jbnd,kbnd,rho0,udir
#ifdef _CUDA
   use m_readinfile, only : ntx,nty,ntz
#endif
   use m_wtime

   use m_boundary_i_inflow_kernel
   use m_boundary_i_inflow_edges

   use m_boundary_i_periodic_kernel
   use m_boundary_j_periodic_kernel
   use m_boundary_k_periodic_kernel

   use m_boundary_j_closed_kernel
   use m_boundary_k_closed_kernel

   use m_boundary_closed_edges

   implicit none
   real, intent(inout):: f1(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout):: f2(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(inout):: tracer(ntracer,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: uvel(nz)
#ifdef _CUDA
   attributes(device) :: f1
   attributes(device) :: f2
   attributes(device) :: tracer
   attributes(device) :: uvel
#endif
   integer, parameter :: icpu=11
#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif

   call cpustart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! i-inflow boundary condition
   if (ibnd==1) then
#ifdef _CUDA
      tx=1;   bx=1
      ty=8;   by=(ny+2+ty-1)/ty
      tz=8;   bz=(nz+2+tz-1)/tz
#endif
      call boundary_i_inflow_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1,uvel,rho0,udir,tracer)

   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Periodic boundary conditions in i-direction.
   if (ibnd==0) then
#ifdef _CUDA
      tx=ntx; bx=(nl  +tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=ntz; bz=(nz+2+tz-1)/tz
#endif
      call boundary_i_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1)
   endif


! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
#ifdef _CUDA
      tx=512; bx=(nl*(nx+2)+tx-1)/tx
      ty=1; by=1
      tz=1; bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1,nl)

      if (ntracer > 0) then
#ifdef _CUDA
      tx=512; bx=(ntracer*(nx+2)+tx-1)/tx
      ty=1; by=1
      tz=1; bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(tracer,ntracer)
      endif
   endif

! Periodic boundary conditions in k-direction.
   if (kbnd==0) then
#ifdef _CUDA
      tx=512; bx=(nl*(nx+2)+tx-1)/tx
      ty=1; by=(ny+2+ty-1)/ty
      tz=1; bz=1
#endif
      call boundary_k_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1)
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Closed boundary conditions in i-direction (at the sides i=1 and i=nx)
   if ((ibnd==11).or.(ibnd==12)) then
      stop 'No-slip bounce-back condition for i=1 not implemented'
   elseif ((ibnd==11).or.(ibnd==21)) then
      stop 'No-slip bounce-back condition for i=nx not implemented'
   elseif ((ibnd==22).or.(ibnd==21)) then
      stop 'Free-slip bounce-back condition for i=1 not implemented'
   elseif ((ibnd==22).or.(ibnd==12)) then
      stop 'Free-slip bounce-back condition for i=nx not implemented'
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! No-slip
   if ((jbnd==11).or.(jbnd==12)) then
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call boundary_j_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,1,-1)
   endif

   if ((jbnd==11).or.(jbnd==21)) then
! No-slip bounce-back condition for j=ny
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call boundary_j_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,ny,-1)
   endif

   if ((kbnd==11).or.(kbnd==12)) then
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
           &(f1,f2,1,-1)
   endif

   if ((kbnd==11).or.(kbnd==21)) then
! No-slip bounce-back condition for k=nz
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=1; bz=1
#endif
      call boundary_k_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,nz,-1)
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Free slip


   if ((jbnd==22).or.(jbnd==21)) then
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call boundary_j_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,1,1)
   endif

   if ((jbnd==22).or.(jbnd==12)) then
! Free-slip bounce-back condition for j=ny
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+tz-1)/tz
#endif
      call boundary_j_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,ny,1)
   endif

   if ((kbnd==22).or.(kbnd==21)) then
! Free-slip bounce-back condition for k=1
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=1; bz=1
#endif
      call boundary_k_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,1,1)
   endif

   if ((kbnd==22).or.(kbnd==12)) then
#ifdef _CUDA
      tx=ntx; bx=(nx+tx-1)/tx
      ty=nty; by=(ny+ty-1)/ty
      tz=1; bz=1
#endif
      call boundary_k_closed_kernel&
#ifdef _CUDA
           &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
           &(f1,f2,nz,1)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update edges for inflow conditions
   if (ibnd==1) then
      call boundary_i_inflow_edges(f1,f2)
   endif

! Update edges for inflow conditions for periodic boundary conditions in i-direction
   if (ibnd==0) then
#ifdef _CUDA
      tx=ntx; bx=(nl  +tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=ntz; bz=(nz+2+tz-1)/tz
#endif
      call boundary_i_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1)
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
#ifdef _CUDA
      tx=512; bx=(nl*(nx+2)+tx-1)/tx
      ty=1; by=1
      tz=1; bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1,nl)

      if (ntracer > 0) then
#ifdef _CUDA
      tx=512; bx=(ntracer*(nx+2)+tx-1)/tx
      ty=1; by=1
      tz=1; bz=(nz+2+tz-1)/tz
#endif
      call boundary_j_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(tracer,ntracer)
      endif
   endif

   if (kbnd==0) then
#ifdef _CUDA
      tx=512; bx=(nl*(nx+2)+tx-1)/tx
      ty=1; by=(ny+2+ty-1)/ty
      tz=1; bz=1
#endif
      call boundary_k_periodic_kernel&
#ifdef _CUDA
      &<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>&
#endif
      &(f1)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bounce-back edges in case of closed boundaries in both y and z directions and when at least one direction has no slip boundaries
   if ((jbnd == 11 .and. kbnd > 10) .or. (jbnd > 10 .and. kbnd == 11)) then
      call boundary_closed_edges(f1,f2,-1)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Freeslip bounce-back edges in case of closed free-slip boundaries in both y and z directions
   if  (jbnd == 22 .and. kbnd == 22) then
      call boundary_closed_edges(f1,f2,1)
   endif

   call cpufinish(icpu)

end subroutine
end module
