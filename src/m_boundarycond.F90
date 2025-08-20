module m_boundarycond
contains
subroutine boundarycond(f,uvel)
   use mod_dimensions
   use mod_D3Q27setup, only : cxs,cys,czs,nl
   use m_readinfile,   only : ibnd,jbnd,kbnd,ntx,nty,ntz
   use m_wtime

   use m_boundary_iinflow
   use m_boundary_iinflow_edges

   use m_boundary_iperiodic
   use m_boundary_jperiodic
   use m_boundary_kperiodic

   use m_boundary_noslipbb_j1_A
   use m_boundary_noslipbb_j1_B
   use m_boundary_noslipbb_jny_A
   use m_boundary_noslipbb_jny_B
   use m_boundary_noslipbb_k1_A
   use m_boundary_noslipbb_k1_B
   use m_boundary_noslipbb_knz_A
   use m_boundary_noslipbb_knz_B
   use m_boundary_noslip_edges

   use m_boundary_freeslipbb_j1_A
   use m_boundary_freeslipbb_j1_B
   use m_boundary_freeslipbb_jny_A
   use m_boundary_freeslipbb_jny_B
   use m_boundary_freeslipbb_k1_A
   use m_boundary_freeslipbb_k1_B
   use m_boundary_freeslipbb_knz_A
   use m_boundary_freeslipbb_knz_B
   use m_boundary_freeslip_edges

   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: uvel(nz)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uvel
#endif
   !integer i,j,k,l,m
   integer, parameter :: icpu=11
#ifdef _CUDA
   integer :: tx, ty, tz, bx, by, bz
#endif

   call cpustart()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (ibnd==1) then
      call boundary_iinflow(f,uvel)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Periodic boundary conditions in i-direction.
   if (ibnd==0) then
      call boundary_iperiodic(f)
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
      call boundary_jperiodic(f)
   endif

! Periodic boundary conditions in k-direction.
   if (kbnd==0) then
      call boundary_kperiodic(f)
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
! Closed boundary conditions in j-direction (at the sides j=1 and j=ny)
   if ((jbnd==11).or.(jbnd==12)) then
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+2+tz-1)/tz
      call boundary_noslipbb_j1_A<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_noslipbb_j1_B<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#else
      call boundary_noslipbb_j1_A                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_noslipbb_j1_B                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#endif
   endif

   if ((jbnd==22).or.(jbnd==21)) then
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+2+tz-1)/tz
      call boundary_freeslipbb_j1_A<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_freeslipbb_j1_B<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#else
      call boundary_freeslipbb_j1_A                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_freeslipbb_j1_B                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#endif
   endif

   if ((jbnd==11).or.(jbnd==21)) then
! No-slip bounce-back condition for j=ny
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+2+tz-1)/tz
      call boundary_noslipbb_jny_A<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_noslipbb_jny_B<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#else
      call boundary_noslipbb_jny_A                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_noslipbb_jny_B                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#endif
   endif

   if ((jbnd==22).or.(jbnd==12)) then
! Free-slip bounce-back condition for j=ny
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=1; by=1
      tz=ntz; bz=(nz+2+tz-1)/tz
      call boundary_freeslipbb_jny_A<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_freeslipbb_jny_B<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#else
      call boundary_freeslipbb_jny_A                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_freeslipbb_jny_B                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#endif
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Closed boundary conditions in k-direction (at the sides k=1 and k=nz)

   if ((kbnd==11).or.(kbnd==12)) then
! No-slip bounce-back condition for k=1
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=1; bz=1
      call boundary_noslipbb_k1_A<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_noslipbb_k1_B<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#else
      call boundary_noslipbb_k1_A                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_noslipbb_k1_B                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#endif
   endif

   if ((kbnd==22).or.(kbnd==21)) then
! Free-slip bounce-back condition for k=1
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=1; bz=1
      call boundary_freeslipbb_k1_A<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_freeslipbb_k1_B<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#else
      call boundary_freeslipbb_k1_A                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_freeslipbb_k1_B                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#endif
   endif

   if ((kbnd==11).or.(kbnd==21)) then
! No-slip bounce-back condition for k=nz
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=1; bz=1
      call boundary_noslipbb_knz_A<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_noslipbb_knz_B<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#else
      call boundary_noslipbb_knz_A                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_noslipbb_knz_B                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#endif
   endif


   if ((kbnd==22).or.(kbnd==12)) then
#ifdef _CUDA
      tx=ntx; bx=(nx+2+tx-1)/tx
      ty=nty; by=(ny+2+ty-1)/ty
      tz=1; bz=1
      call boundary_freeslipbb_knz_A<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_freeslipbb_knz_B<<<dim3(bx,by,bz), dim3(tx,ty,tz)>>>(f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#else
      call boundary_freeslipbb_knz_A                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
      call boundary_freeslipbb_knz_B                                    (f,nx+2,ny+2,nz+2,nl,cxs,cys,czs)
#endif
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update edges for inflow conditions
   if (ibnd==1) then
      call boundary_iinflow_edges(f)
   endif

! Update edges for inflow conditions for periodic boundary conditions in i-direction
   if (ibnd==0) then
      call boundary_iperiodic(f)
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
      call boundary_jperiodic(f)
   endif

   if (kbnd==0) then
      call boundary_kperiodic(f)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bounce-back edges in case of closed boundaries in both y and z directions and when at least one direction has no slip boundaries
   if ((jbnd == 11 .and. kbnd > 10) .or. (jbnd > 10 .and. kbnd == 11)) then
      call boundary_noslip_edges(f)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Freeslip bounce-back edges in case of closed free-slip boundaries in both y and z directions
   if  (jbnd == 22 .and. kbnd == 22) then
      call boundary_freeslip_edges(f)
   endif

   call cpufinish(icpu)

end subroutine
end module
