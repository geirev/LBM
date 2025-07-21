module m_boundarycond
contains
subroutine boundarycond(f,rho,u,v,w,uvel)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_wtime

   use m_boundary_noslipbb_j1_A
   use m_boundary_noslipbb_j1_B
   use m_boundary_noslipbb_jny_A
   use m_boundary_noslipbb_jny_B
   use m_boundary_noslipbb_k1_A
   use m_boundary_noslipbb_k1_B
   use m_boundary_noslipbb_knz_A
   use m_boundary_noslipbb_knz_B

   use m_boundary_freeslipbb_j1_A
   use m_boundary_freeslipbb_j1_B
   use m_boundary_freeslipbb_jny_A
   use m_boundary_freeslipbb_jny_B
   use m_boundary_freeslipbb_k1_A
   use m_boundary_freeslipbb_k1_B
   use m_boundary_freeslipbb_knz_A
   use m_boundary_freeslipbb_knz_B

   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: rho(nx,ny,nz)
   real, intent(in)   :: u(nx,ny,nz)
   real, intent(in)   :: v(nx,ny,nz)
   real, intent(in)   :: w(nx,ny,nz)
   real, intent(in)   :: uvel(nz)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: uvel
#endif
   real tmp
   real tmpval
   integer i,j,k,l,m,ja,ka,ip,ix,jy,kz
   integer, parameter :: icpu=11
   real, parameter :: pi=3.1415927410125732
   integer :: tx, ty, tz, bx, by, bz


   call cpustart()


! Periodic boundary conditions in i-direction.
   if (ibnd==0) then
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
      do l=1,nl
         f(l,0,j,k)   =f(l,nx,j,k)
         f(l,nx+1,j,k)=f(l,1,j,k)
      enddo
      enddo
      enddo
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do i=0,nx+1
      do l=1,nl
         f(l,i,0,k)   =f(l,i,ny,k)
         f(l,i,ny+1,k)=f(l,i,1,k)
      enddo
      enddo
      enddo
   endif

! Periodic boundary conditions in k-direction.
   if (kbnd==0) then
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do j=0,ny+1
      do i=0,nx+1
      do l=1,nl
         f(l,i,j,0)   =f(l,i,j,nz)
         f(l,i,j,nz+1)=f(l,i,j,1)
      enddo
      enddo
      enddo
   endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inflow outflow boundary conditions in i-direction
   if (ibnd==1) then
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>> 
#endif
      do k=0,nz+1
      do j=0,ny+1
         ka = min(max(k,1), nz)
         ja = min(max(j,1), ny)
         do l=1,nl
            f(l,0,j,k) = f(l,1,j,k) - 2.0*weights(l)*rho0*(cxs(l)*uvel(ka)*cos(udir*pi/180.0)&
                                                         + cys(l)*uvel(ka)*sin(udir*pi/180.0)&
                                                         + czs(l)*0.0)/cs2
         enddo
         do l=2,nl-1,2
            tmp=f(l,0,j,k)
            if (cxs(l)==1)   f(l,0,j,k)=f(l+1,0,j,k)
            if (cxs(l+1)==1) f(l+1,0,j,k)=tmp
         enddo
      enddo
      enddo


#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
         f(:,nx+1,j,k)=f(:,nx,j,k)
      enddo
      enddo

   elseif ((ibnd==11).or.(ibnd==12)) then
! No-slip bounce-back condition for i=1
      stop
   elseif ((ibnd==11).or.(ibnd==21)) then
! No-slip bounce-back condition for i=nx
      stop
   elseif ((ibnd==22).or.(ibnd==21)) then
! Free-slip bounce-back condition for i=1
      stop
   elseif ((ibnd==22).or.(ibnd==12)) then
! Free-slip bounce-back condition for i=nx
      stop
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Closed boundary conditions in j-direction
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
! Closed boundary conditions in k-direction

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

!                           1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7
!   integer :: cxs(1:nl) = [0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1]
!   integer :: cys(1:nl) = [0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 0, 0,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1]
!   integer :: czs(1:nl) = [0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1]


! Periodic boundary conditions in i-direction.
   if (ibnd==0) then
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
      do l=1,nl
         f(l,0,j,k)   =f(l,nx,j,k)
         f(l,nx+1,j,k)=f(l,1,j,k)
      enddo
      enddo
      enddo
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do k=0,nz+1
      do i=0,nx+1
      do l=1,nl
         f(l,i,0,k)   =f(l,i,ny,k)
         f(l,i,ny+1,k)=f(l,i,1,k)
      enddo
      enddo
      enddo
   endif

! Periodic boundary conditions in k-direction.
   if (kbnd==0) then
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
      do j=0,ny+1
      do i=0,nx+1
      do l=1,nl
         f(l,i,j,0)   =f(l,i,j,nz)
         f(l,i,j,nz+1)=f(l,i,j,1)
      enddo
      enddo
      enddo
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bounce-back edges in case of closed boundaries in both y and z directions
! Only needed in case of closed boundaries in both y and z directions

! The free-slip edges have not yet been implemented and I didn't care about the corners yet.
   if  (jbnd > 10 .and. kbnd > 10) then

! === Edge along j=1 and k=1 === No slip
#ifdef CUDAS
!$cuf kernel do(1) <<<*,*>>>
#endif
    do i=1,nx
       do l=1,nl
          ip = bounce(l)
          if (cys(l) < 0.0 .and. czs(l) < 0.0) then
             tmp=f(l,i,0,0); f(ip,i,0,0)=tmp
             ix=i-cxs(l); jy=0-cys(l); kz=0-czs(l); tmp=f(l, ix, jy, kz); f(l, i, 0, 0)=tmp
          endif
          if (cys(l) <= 0.0 .and. czs(l) < 0.0) then
             tmp=f(l,i,1,0); f(ip,i,1,0) = tmp
             ix=i-cxs(l); jy=1-cys(l); kz=0-czs(l); tmp=f(l, ix, jy, kz); f(l, i, 1, 0)=tmp
          endif
          if (cys(l) < 0.0 .and. czs(l) <= 0.0) then
             tmp=f(l,i,0,1); f(ip,i,0,1) = tmp
             ix=i-cxs(l); jy=0-cys(l); kz=1-czs(l); tmp=f(l, ix, jy, kz); f(l, i, 0, 1)=tmp
          endif
       enddo
    enddo

! === Edge along j=ny and k=1 === No slip
#ifdef CUDAS
!$cuf kernel do(1) <<<*,*>>>
#endif
    do i=1,nx
       do l=1,nl
          ip = bounce(l)
          if (cys(l) > 0.0 .and. czs(l) < 0.0) then
            tmp=f(l,i,ny+1,0);                               f(ip,i,ny+1,0)=tmp
            ix=i-cxs(l); jy=ny+1-cys(l); kz=0-czs(l); tmp=f(l, ix, jy, kz); f(l, i, ny+1, 0)=tmp
          endif
          if (cys(l) >= 0.0 .and. czs(l) > 0.0) then
             tmp=f(l,i,ny,0);                               f(ip,i,ny,0)=tmp
             ix=i-cxs(l); jy=ny-cys(l); kz=0-czs(l); tmp=f(l, ix, jy, kz); f(l, i, ny, 0)=tmp
          endif
          if (cys(l) > 0.0 .and. czs(l) >= 0.0) then
             tmp=f(l,i,ny+1,1);                              f(ip,i,ny+1,1)=tmp
             ix=i-cxs(l); jy=ny-cys(l); kz=0-czs(l); tmp=f(l,ix,jy,kz); f( l,i,ny+1,1)=tmp
          endif
       enddo
    enddo

! === Edge along j=1 and k=nz === No slip
#ifdef CUDAS
!$cuf kernel do(1) <<<*,*>>>
#endif
       do i=1,nx
          do l=1,nl
             ip = bounce(l)
             if (cys(l) < 0.0 .and. czs(l) > 0.0) then
               tmp=f(l,i,0,nz+1);                               f(ip,i,0,nz+1)=tmp
               ix=i-cxs(l); jy=0-cys(l); kz=nz+1-czs(l); tmp=f(l,ix,jy,kz); f(l, i,0,nz+1)=tmp
             endif
             if (cys(l) <= 0.0 .and. czs(l) > 0.0) then
                tmp=f(l,i,1,nz+1);                              f(ip,i,1,nz+1)=tmp
                ix=i-cxs(l); jy=1-cys(l); kz=nz+1-czs(l); tmp=f(l,ix,jy,kz); f( l,i,1,nz+1)=tmp
             endif
             if (cys(l) < 0.0 .and. czs(l) >= 0.0) then
                tmp=f(l,i,0,nz);                              f(ip,i,0,nz)=tmp
                ix=i-cxs(l); jy=0-cys(l); kz=nz-czs(l); tmp=f(l,ix,jy,kz); f( l,i,0,nz)=tmp
             endif
          enddo
       enddo

! === Edge along j=ny and k=nz === No slip
#ifdef CUDAS
!$cuf kernel do(1) <<<*,*>>>
#endif
       do i=1,nx
          do l=1,nl
             ip = bounce(l)
             if (cys(l) > 0.0 .and. czs(l) > 0.0) then
               tmp=f(l,i,ny+1,nz+1);                               f(ip,i,ny+1,nz+1)=tmp
               ix=i-cxs(l); jy=ny+1-cys(l); kz=nz+1-czs(l); tmp=f(l,ix,jy,kz); f(l, i,ny+1,nz+1)=tmp
             endif
             if (cys(l) >= 0.0 .and. czs(l) > 0.0) then
                tmp=f(l,i,ny,nz+1);                               f(ip,i,ny,nz+1)=tmp
                ix=i-cxs(l); jy=ny-cys(l); kz=nz+1-czs(l); tmp=f(l,ix,jy,kz); f( l,i,ny,nz+1)=tmp
             endif
             if (cys(l) > 0.0 .and. czs(l) >= 0.0) then
                tmp=f(l,i,ny+1,nz);                              f(ip,i,ny+1,nz)=tmp
                ix=i-cxs(l); jy=ny+1-cys(l); kz=nz-czs(l); tmp=f(l,ix,jy,kz); f( l,i,ny+1,nz)=tmp
             endif
          enddo
       enddo

    endif

   call cpufinish(icpu)

end subroutine
end module
