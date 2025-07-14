module m_boundarycond
contains
subroutine boundarycond(f,rho,u,v,w,uvel)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_wtime
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


   call cpustart()

! Periodic boundary conditions in i-direction.
   if (ibnd==0) then
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
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
!$cuf kernel do(2) <<<*,*>>>
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
!$cuf kernel do(2) <<<*,*>>>
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
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
         ka = min(max(k,1), nz)
         ja = min(max(j,1), ny)
!         ka=k; if (k<1) ka=1; if (k>nz) ka = nz
!         ja=j; if (j<1) ja=1; if (j>ny) ja = ny
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
!$cuf kernel do(2) <<<*,*>>>
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
! No-slip bounce-back condition for j=1
!     print *,'No-slip bounce-back condition for j=1'
!     do l = 1, nl
!        if (cys(l) < 0.0) then
!           do m = 1, nl
!              if (cxs(m) == -cxs(l) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
!                f(m, i, j, nz+1) = f(l, i, j, nz)
!                print '(5(a,i2),a)','f(',m,',i,0,k) = f(',l,',i-cxs(',l,'),1,k-czs(',l,'))'
!                exit
!              endif
!           enddo
!        endif
!     enddo
#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=1,nz
      do i=1,nx
         tmpval=f( 5,i,0,k); f( 4,i,0,k)=tmpval
         tmpval=f( 9,i,0,k); f( 8,i,0,k)=tmpval
         tmpval=f(10,i,0,k); f(11,i,0,k)=tmpval
         tmpval=f(15,i,0,k); f(14,i,0,k)=tmpval
         tmpval=f(18,i,0,k); f(19,i,0,k)=tmpval
         tmpval=f(21,i,0,k); f(20,i,0,k)=tmpval
         tmpval=f(22,i,0,k); f(23,i,0,k)=tmpval
         tmpval=f(25,i,0,k); f(24,i,0,k)=tmpval
         tmpval=f(27,i,0,k); f(26,i,0,k)=tmpval

         ix=i-cxs( 5); kz=k-czs( 5); tmpval=f( 5,ix,1,kz); f( 5,i,0,k)=tmpval
         ix=i-cxs( 9); kz=k-czs( 9); tmpval=f( 9,ix,1,kz); f( 9,i,0,k)=tmpval
         ix=i-cxs(10); kz=k-czs(10); tmpval=f(10,ix,1,kz); f(10,i,0,k)=tmpval
         ix=i-cxs(15); kz=k-czs(15); tmpval=f(15,ix,1,kz); f(15,i,0,k)=tmpval
         ix=i-cxs(18); kz=k-czs(18); tmpval=f(18,ix,1,kz); f(18,i,0,k)=tmpval
         ix=i-cxs(21); kz=k-czs(21); tmpval=f(21,ix,1,kz); f(21,i,0,k)=tmpval
         ix=i-cxs(22); kz=k-czs(22); tmpval=f(22,ix,1,kz); f(22,i,0,k)=tmpval
         ix=i-cxs(25); kz=k-czs(25); tmpval=f(25,ix,1,kz); f(25,i,0,k)=tmpval
         ix=i-cxs(27); kz=k-czs(27); tmpval=f(27,ix,1,kz); f(27,i,0,k)=tmpval
      enddo
      enddo
!      do k=1,nz+1
!      do i=1,nx+1
!         f( 4,i,0,k) = f( 5,i-cxs( 5),1,k-czs( 5))
!         f( 8,i,0,k) = f( 9,i-cxs( 9),1,k-czs( 9))
!         f(11,i,0,k) = f(10,i-cxs(10),1,k-czs(10))
!         f(14,i,0,k) = f(15,i-cxs(15),1,k-czs(15))
!         f(19,i,0,k) = f(18,i-cxs(18),1,k-czs(18))
!         f(20,i,0,k) = f(21,i-cxs(21),1,k-czs(21))
!         f(23,i,0,k) = f(22,i-cxs(22),1,k-czs(22))
!         f(24,i,0,k) = f(25,i-cxs(25),1,k-czs(25))
!         f(26,i,0,k) = f(27,i-cxs(27),1,k-czs(27))
!      enddo
!      enddo
   endif

   if ((jbnd==22).or.(jbnd==21)) then
! Free-slip bounce-back condition for j=1
!     do l = 1, nl
!        if (cys(l) < 0.0) then
!           do m = 1, nl
!              if (cxs(m) ==  cxs(l) .and. cys(m) ==  -cys(l) .and. czs(m) == czs(l)) then
!                print '(5(a,i2),a)','f(',m,',i,0,k) = f(',l,',i-cxs(',l,'),1,k-czs(',l,'))'
!                exit
!              endif
!           enddo
!        endif
!     enddo
#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=1,nz
      do i=1,nx
         tmpval=f( 5,i,0,k); f( 4,i,0,k)=tmpval
         tmpval=f( 9,i,0,k); f(11,i,0,k)=tmpval
         tmpval=f(10,i,0,k); f( 8,i,0,k)=tmpval
         tmpval=f(15,i,0,k); f(19,i,0,k)=tmpval
         tmpval=f(18,i,0,k); f(14,i,0,k)=tmpval
         tmpval=f(21,i,0,k); f(24,i,0,k)=tmpval
         tmpval=f(22,i,0,k); f(26,i,0,k)=tmpval
         tmpval=f(25,i,0,k); f(20,i,0,k)=tmpval
         tmpval=f(27,i,0,k); f(23,i,0,k)=tmpval

         ix=i-cxs( 5); kz=k-czs( 5); tmpval=f( 5,ix,1,kz); f( 5,i,0,k)=tmpval
         ix=i-cxs( 9); kz=k-czs( 9); tmpval=f( 9,ix,1,kz); f( 9,i,0,k)=tmpval
         ix=i-cxs(10); kz=k-czs(10); tmpval=f(10,ix,1,kz); f(10,i,0,k)=tmpval
         ix=i-cxs(15); kz=k-czs(15); tmpval=f(15,ix,1,kz); f(15,i,0,k)=tmpval
         ix=i-cxs(18); kz=k-czs(18); tmpval=f(18,ix,1,kz); f(18,i,0,k)=tmpval
         ix=i-cxs(21); kz=k-czs(21); tmpval=f(21,ix,1,kz); f(21,i,0,k)=tmpval
         ix=i-cxs(22); kz=k-czs(22); tmpval=f(22,ix,1,kz); f(22,i,0,k)=tmpval
         ix=i-cxs(25); kz=k-czs(25); tmpval=f(25,ix,1,kz); f(25,i,0,k)=tmpval
         ix=i-cxs(27); kz=k-czs(27); tmpval=f(27,ix,1,kz); f(27,i,0,k)=tmpval
      enddo
      enddo
!      do k=1,nz
!      do i=1,nx
!         f( 4,i,0,k) = f( 5,i-cxs( 5),1,k-czs( 5))
!         f(11,i,0,k) = f( 9,i-cxs( 9),1,k-czs( 9))
!         f( 8,i,0,k) = f(10,i-cxs(10),1,k-czs(10))
!         f(19,i,0,k) = f(15,i-cxs(15),1,k-czs(15))
!         f(14,i,0,k) = f(18,i-cxs(18),1,k-czs(18))
!         f(24,i,0,k) = f(21,i-cxs(21),1,k-czs(21))
!         f(26,i,0,k) = f(22,i-cxs(22),1,k-czs(22))
!         f(20,i,0,k) = f(25,i-cxs(25),1,k-czs(25))
!         f(23,i,0,k) = f(27,i-cxs(27),1,k-czs(27))
!      enddo
!      enddo

   endif

   if ((jbnd==11).or.(jbnd==21)) then
! No-slip bounce-back condition for j=ny
!     print *,'No-slip bounce-back condition for j=ny'
!     do l = 1, nl
!        if (cys(l) > 0.0) then
!           do m = 1, nl
!              if (cxs(m) ==  -cxs(l) .and. cys(m) ==  -cys(l) .and. czs(m) == -czs(l)) then
!                print '(5(a,i2),a)','f(',m,',i,ny+1,k) = f(',l,',i-cxs(',l,'),ny,k-czs(',l,'))'
!                exit
!              endif
!           enddo
!        endif
!     enddo
!     stop
#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=1,nz
      do i=1,nx
         tmpval=f( 4,i,ny+1,k); f( 5,i,ny+1,k)=tmpval
         tmpval=f( 8,i,ny+1,k); f( 9,i,ny+1,k)=tmpval
         tmpval=f(11,i,ny+1,k); f(10,i,ny+1,k)=tmpval
         tmpval=f(14,i,ny+1,k); f(15,i,ny+1,k)=tmpval
         tmpval=f(19,i,ny+1,k); f(18,i,ny+1,k)=tmpval
         tmpval=f(20,i,ny+1,k); f(21,i,ny+1,k)=tmpval
         tmpval=f(23,i,ny+1,k); f(22,i,ny+1,k)=tmpval
         tmpval=f(24,i,ny+1,k); f(25,i,ny+1,k)=tmpval
         tmpval=f(26,i,ny+1,k); f(27,i,ny+1,k)=tmpval

         ix=i-cxs( 4); kz=k-czs( 4); tmpval=f( 4,ix,ny,kz); f( 4,i,ny+1,k)=tmpval
         ix=i-cxs( 8); kz=k-czs( 8); tmpval=f( 8,ix,ny,kz); f( 8,i,ny+1,k)=tmpval
         ix=i-cxs(11); kz=k-czs(11); tmpval=f(11,ix,ny,kz); f(11,i,ny+1,k)=tmpval
         ix=i-cxs(14); kz=k-czs(14); tmpval=f(14,ix,ny,kz); f(14,i,ny+1,k)=tmpval
         ix=i-cxs(19); kz=k-czs(19); tmpval=f(19,ix,ny,kz); f(19,i,ny+1,k)=tmpval
         ix=i-cxs(20); kz=k-czs(20); tmpval=f(20,ix,ny,kz); f(20,i,ny+1,k)=tmpval
         ix=i-cxs(23); kz=k-czs(23); tmpval=f(23,ix,ny,kz); f(23,i,ny+1,k)=tmpval
         ix=i-cxs(24); kz=k-czs(24); tmpval=f(24,ix,ny,kz); f(24,i,ny+1,k)=tmpval
         ix=i-cxs(26); kz=k-czs(26); tmpval=f(26,ix,ny,kz); f(26,i,ny+1,k)=tmpval
      enddo
      enddo
!      do k=1,nz
!      do i=1,nx
!         f( 5,i,ny+1,k) = f( 4,i-cxs( 4),ny,k-czs( 4))
!         f( 9,i,ny+1,k) = f( 8,i-cxs( 8),ny,k-czs( 8))
!         f(10,i,ny+1,k) = f(11,i-cxs(11),ny,k-czs(11))
!         f(15,i,ny+1,k) = f(14,i-cxs(14),ny,k-czs(14))
!         f(18,i,ny+1,k) = f(19,i-cxs(19),ny,k-czs(19))
!         f(21,i,ny+1,k) = f(20,i-cxs(20),ny,k-czs(20))
!         f(22,i,ny+1,k) = f(23,i-cxs(23),ny,k-czs(23))
!         f(25,i,ny+1,k) = f(24,i-cxs(24),ny,k-czs(24))
!         f(27,i,ny+1,k) = f(26,i-cxs(26),ny,k-czs(26))
!      enddo
!      enddo
   endif

   if ((jbnd==22).or.(jbnd==12)) then
! Free-slip bounce-back condition for j=ny
!     do l = 1, nl
!        if (cys(l) > 0.0) then
!           do m = 1, nl
!              if (cxs(m) == cxs(l) .and. cys(m) == -cys(l) .and. czs(m) == czs(l)) then
!                print '(5(a,i2),a)','f(',m,',i,ny+1,k) = f(',l,',i-cxs(',l,'),ny,k-czs(',l,'))'
!                exit
!              endif
!           enddo
!        endif
!     enddo
#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=1,nz
      do i=1,nx
         tmpval=f( 4,i,ny+1,k); f( 5,i,ny+1,k)=tmpval
         tmpval=f( 8,i,ny+1,k); f(10,i,ny+1,k)=tmpval
         tmpval=f(11,i,ny+1,k); f( 9,i,ny+1,k)=tmpval
         tmpval=f(14,i,ny+1,k); f(18,i,ny+1,k)=tmpval
         tmpval=f(19,i,ny+1,k); f(15,i,ny+1,k)=tmpval
         tmpval=f(20,i,ny+1,k); f(25,i,ny+1,k)=tmpval
         tmpval=f(23,i,ny+1,k); f(27,i,ny+1,k)=tmpval
         tmpval=f(24,i,ny+1,k); f(21,i,ny+1,k)=tmpval
         tmpval=f(26,i,ny+1,k); f(22,i,ny+1,k)=tmpval

         ix=i-cxs( 4); kz=k-czs( 4); tmpval=f( 4,ix,ny,kz); f( 4,i,ny+1,k)=tmpval
         ix=i-cxs( 8); kz=k-czs( 8); tmpval=f( 8,ix,ny,kz); f( 8,i,ny+1,k)=tmpval
         ix=i-cxs(11); kz=k-czs(11); tmpval=f(11,ix,ny,kz); f(11,i,ny+1,k)=tmpval
         ix=i-cxs(14); kz=k-czs(14); tmpval=f(14,ix,ny,kz); f(14,i,ny+1,k)=tmpval
         ix=i-cxs(19); kz=k-czs(19); tmpval=f(19,ix,ny,kz); f(19,i,ny+1,k)=tmpval
         ix=i-cxs(20); kz=k-czs(20); tmpval=f(20,ix,ny,kz); f(20,i,ny+1,k)=tmpval
         ix=i-cxs(23); kz=k-czs(23); tmpval=f(23,ix,ny,kz); f(23,i,ny+1,k)=tmpval
         ix=i-cxs(24); kz=k-czs(24); tmpval=f(24,ix,ny,kz); f(24,i,ny+1,k)=tmpval
         ix=i-cxs(26); kz=k-czs(26); tmpval=f(26,ix,ny,kz); f(26,i,ny+1,k)=tmpval
      enddo
      enddo
!      do k=1,nz
!      do i=1,nx
!         f( 5,i,ny+1,k) = f( 4,i-cxs( 4),ny,k-czs( 4))
!         f(10,i,ny+1,k) = f( 8,i-cxs( 8),ny,k-czs( 8))
!         f( 9,i,ny+1,k) = f(11,i-cxs(11),ny,k-czs(11))
!         f(18,i,ny+1,k) = f(14,i-cxs(14),ny,k-czs(14))
!         f(15,i,ny+1,k) = f(19,i-cxs(19),ny,k-czs(19))
!         f(25,i,ny+1,k) = f(20,i-cxs(20),ny,k-czs(20))
!         f(27,i,ny+1,k) = f(23,i-cxs(23),ny,k-czs(23))
!         f(21,i,ny+1,k) = f(24,i-cxs(24),ny,k-czs(24))
!         f(22,i,ny+1,k) = f(26,i-cxs(26),ny,k-czs(26))
!      enddo
!      enddo
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Closed boundary conditions in k-direction


   if ((kbnd==11).or.(kbnd==12)) then
! No-slip bounce-back condition for k=1
! Compute drift into ghost nodes followed by switch to opposite direction
!     do l = 1, nl
!        if (czs(l) < 0.0) then
!           do m = 1, nl
!              if (cxs(m) == -cxs(l) .and. cys(m) == -cys(l) .and. czs(m) == -czs(l)) then
!                f(m, i, j, nz+1) = f(l, i, j, nz)
!                print '(5(a,i2),a)','f(',m,',i,j,0) = f(',l,',i-cxs(',l,'),j-cys(',l,'),1)'
!                exit
!              endif
!           enddo
!        endif
!     enddo
#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
#endif
      do j=1,ny
      do i=1,nx
         tmpval=f( 6,i,j,0); f( 7,i,j,0)=tmpval
         tmpval=f(12,i,j,0); f(13,i,j,0)=tmpval
         tmpval=f(15,i,j,0); f(14,i,j,0)=tmpval
         tmpval=f(17,i,j,0); f(16,i,j,0)=tmpval
         tmpval=f(19,i,j,0); f(18,i,j,0)=tmpval
         tmpval=f(21,i,j,0); f(20,i,j,0)=tmpval
         tmpval=f(22,i,j,0); f(23,i,j,0)=tmpval
         tmpval=f(24,i,j,0); f(25,i,j,0)=tmpval
         tmpval=f(26,i,j,0); f(27,i,j,0)=tmpval

         ix=i-cxs( 6); jy=j-cys( 6); tmpval=f( 6,ix,jy,1); f( 6,i,jy+1,0)=tmpval
         ix=i-cxs(12); jy=j-cys(12); tmpval=f(12,ix,jy,1); f(12,i,jy+1,0)=tmpval
         ix=i-cxs(15); jy=j-cys(15); tmpval=f(15,ix,jy,1); f(15,i,jy+1,0)=tmpval
         ix=i-cxs(17); jy=j-cys(17); tmpval=f(17,ix,jy,1); f(17,i,jy+1,0)=tmpval
         ix=i-cxs(19); jy=j-cys(19); tmpval=f(19,ix,jy,1); f(19,i,jy+1,0)=tmpval
         ix=i-cxs(21); jy=j-cys(21); tmpval=f(21,ix,jy,1); f(21,i,jy+1,0)=tmpval
         ix=i-cxs(22); jy=j-cys(22); tmpval=f(22,ix,jy,1); f(22,i,jy+1,0)=tmpval
         ix=i-cxs(24); jy=j-cys(24); tmpval=f(24,ix,jy,1); f(24,i,jy+1,0)=tmpval
         ix=i-cxs(26); jy=j-cys(26); tmpval=f(26,ix,jy,1); f(26,i,jy+1,0)=tmpval
      enddo
      enddo
!      do j=0,ny
!      do i=0,nx
!         f( 7,i,j,0)=f( 6,i-cxs( 6),j-cys( 6),1)
!         f(13,i,j,0)=f(12,i-cxs(12),j-cys(12),1)
!         f(14,i,j,0)=f(15,i-cxs(15),j-cys(15),1)
!         f(16,i,j,0)=f(17,i-cxs(17),j-cys(17),1)
!         f(18,i,j,0)=f(19,i-cxs(19),j-cys(19),1)
!         f(20,i,j,0)=f(21,i-cxs(21),j-cys(21),1)
!         f(23,i,j,0)=f(22,i-cxs(22),j-cys(22),1)
!         f(25,i,j,0)=f(24,i-cxs(24),j-cys(24),1)
!         f(27,i,j,0)=f(26,i-cxs(26),j-cys(26),1)
!      enddo
!      enddo
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if ((kbnd==22).or.(kbnd==21)) then
! Free-slip bounce-back condition for k=1
!         do l = 1, nl
!            if (czs(l) < 0.0) then
!               do m = 1, nl
!                  if (cxs(m) ==  cxs(l) .and. cys(m) ==  cys(l) .and. czs(m) == -czs(l)) then
!                    print '(5(a,i2),a)','f(',m,',i,j,0) = f(',l,',i-cxs(',l,'),j-cys(',l,'),1)'
!                    exit
!                  endif
!               enddo
!            endif
!         enddo
!      stop
#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
#endif
      do j=1,ny
      do i=1,nx
         tmpval=f( 6,i,j,0); f( 7,i,j,0)=tmpval
         tmpval=f(12,i,j,0); f(16,i,j,0)=tmpval
         tmpval=f(15,i,j,0); f(18,i,j,0)=tmpval
         tmpval=f(17,i,j,0); f(13,i,j,0)=tmpval
         tmpval=f(19,i,j,0); f(14,i,j,0)=tmpval
         tmpval=f(21,i,j,0); f(27,i,j,0)=tmpval
         tmpval=f(22,i,j,0); f(25,i,j,0)=tmpval
         tmpval=f(24,i,j,0); f(23,i,j,0)=tmpval
         tmpval=f(26,i,j,0); f(20,i,j,0)=tmpval

         ix=i-cxs( 6); jy=j-cys( 6); tmpval=f( 6,ix,jy,1); f( 6,i,jy+1,0)=tmpval
         ix=i-cxs(12); jy=j-cys(12); tmpval=f(12,ix,jy,1); f(12,i,jy+1,0)=tmpval
         ix=i-cxs(15); jy=j-cys(15); tmpval=f(15,ix,jy,1); f(15,i,jy+1,0)=tmpval
         ix=i-cxs(17); jy=j-cys(17); tmpval=f(17,ix,jy,1); f(17,i,jy+1,0)=tmpval
         ix=i-cxs(19); jy=j-cys(19); tmpval=f(19,ix,jy,1); f(19,i,jy+1,0)=tmpval
         ix=i-cxs(21); jy=j-cys(21); tmpval=f(21,ix,jy,1); f(21,i,jy+1,0)=tmpval
         ix=i-cxs(22); jy=j-cys(22); tmpval=f(22,ix,jy,1); f(22,i,jy+1,0)=tmpval
         ix=i-cxs(24); jy=j-cys(24); tmpval=f(24,ix,jy,1); f(24,i,jy+1,0)=tmpval
         ix=i-cxs(26); jy=j-cys(26); tmpval=f(26,ix,jy,1); f(26,i,jy+1,0)=tmpval
      enddo
      enddo
!      do j=1,ny+1
!      do i=1,nx+1
!         f( 7,i,j,0) = f( 6,i-cxs( 6),j-cys( 6),1)
!         f(16,i,j,0) = f(12,i-cxs(12),j-cys(12),1)
!         f(18,i,j,0) = f(15,i-cxs(15),j-cys(15),1)
!         f(13,i,j,0) = f(17,i-cxs(17),j-cys(17),1)
!         f(14,i,j,0) = f(19,i-cxs(19),j-cys(19),1)
!         f(27,i,j,0) = f(21,i-cxs(21),j-cys(21),1)
!         f(25,i,j,0) = f(22,i-cxs(22),j-cys(22),1)
!         f(23,i,j,0) = f(24,i-cxs(24),j-cys(24),1)
!         f(20,i,j,0) = f(26,i-cxs(26),j-cys(26),1)
!      enddo
!      enddo
   endif

   if ((kbnd==11).or.(kbnd==21)) then
! No-slip bounce-back condition for k=nz
!         do l = 1, nl
!            if (czs(l) > 0.0) then
!               do m = 1, nl
!                  if (cxs(m) ==  -cxs(l) .and. cys(m) ==  -cys(l) .and. czs(m) == -czs(l)) then
!                    print '(5(a,i2),a)','f(',m,',i,j,nz+1) = f(',l,',i-cxs(',l,'),j-cys(',l,'),nz)'
!                    exit
!                  endif
!               enddo
!            endif
!         enddo
!         stop
#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
#endif
      do j=1,ny+1
      do i=1,nx+1
         tmpval=f( 7,i,j,nz+1); f( 6,i,j,nz+1)=tmpval
         tmpval=f(13,i,j,nz+1); f(12,i,j,nz+1)=tmpval
         tmpval=f(14,i,j,nz+1); f(15,i,j,nz+1)=tmpval
         tmpval=f(16,i,j,nz+1); f(17,i,j,nz+1)=tmpval
         tmpval=f(18,i,j,nz+1); f(19,i,j,nz+1)=tmpval
         tmpval=f(20,i,j,nz+1); f(21,i,j,nz+1)=tmpval
         tmpval=f(23,i,j,nz+1); f(22,i,j,nz+1)=tmpval
         tmpval=f(25,i,j,nz+1); f(24,i,j,nz+1)=tmpval
         tmpval=f(27,i,j,nz+1); f(26,i,j,nz+1)=tmpval

         ix=i-cxs( 7); jy=j-cys( 7); tmpval=f( 7,ix,jy,nz); f( 7,i,jy+1,nz+1)=tmpval
         ix=i-cxs(13); jy=j-cys(13); tmpval=f(13,ix,jy,nz); f(13,i,jy+1,nz+1)=tmpval
         ix=i-cxs(14); jy=j-cys(14); tmpval=f(14,ix,jy,nz); f(14,i,jy+1,nz+1)=tmpval
         ix=i-cxs(16); jy=j-cys(16); tmpval=f(16,ix,jy,nz); f(16,i,jy+1,nz+1)=tmpval
         ix=i-cxs(18); jy=j-cys(18); tmpval=f(18,ix,jy,nz); f(18,i,jy+1,nz+1)=tmpval
         ix=i-cxs(20); jy=j-cys(20); tmpval=f(20,ix,jy,nz); f(20,i,jy+1,nz+1)=tmpval
         ix=i-cxs(23); jy=j-cys(23); tmpval=f(23,ix,jy,nz); f(23,i,jy+1,nz+1)=tmpval
         ix=i-cxs(25); jy=j-cys(25); tmpval=f(25,ix,jy,nz); f(25,i,jy+1,nz+1)=tmpval
         ix=i-cxs(27); jy=j-cys(27); tmpval=f(27,ix,jy,nz); f(27,i,jy+1,nz+1)=tmpval
      enddo
      enddo

!      do j=1,ny
!      do i=1,nx
!         f( 6,i,j,nz+1) = f( 7,i-cxs( 7),j-cys( 7),nz)
!         f(12,i,j,nz+1) = f(13,i-cxs(13),j-cys(13),nz)
!         f(15,i,j,nz+1) = f(14,i-cxs(14),j-cys(14),nz)
!         f(17,i,j,nz+1) = f(16,i-cxs(16),j-cys(16),nz)
!         f(19,i,j,nz+1) = f(18,i-cxs(18),j-cys(18),nz)
!         f(21,i,j,nz+1) = f(20,i-cxs(20),j-cys(20),nz)
!         f(22,i,j,nz+1) = f(23,i-cxs(23),j-cys(23),nz)
!         f(24,i,j,nz+1) = f(25,i-cxs(25),j-cys(25),nz)
!         f(26,i,j,nz+1) = f(27,i-cxs(27),j-cys(27),nz)
!      enddo
!      enddo
   endif


   if ((kbnd==22).or.(kbnd==12)) then
! Free-slip bounce-back condition for k=nz
!         f(l,i,j,0) = feq(l,i-cxs(l),j-cys(l),nz+1-czs(l))
!         do l = 1, nl
!            if (czs(l) > 0.0) then
!               do m = 1, nl
!                  if (cxs(m) == cxs(l) .and. cys(m) == cys(l) .and. czs(m) == -czs(l)) then
!                    print '(5(a,i2),a)','f(',m,',i,j,nz+1) = f(',l,',i-cxs(',l,'),j-cys(',l,'),nz)'
!                    exit
!                  endif
!               enddo
!            endif
!         enddo
!      stop

#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
#endif
      do j = 1, ny
      do i = 1, nx
         tmpval=f( 7,i,j,nz+1); f( 6,i,j,nz+1)=tmpval
         tmpval=f(13,i,j,nz+1); f(17,i,j,nz+1)=tmpval
         tmpval=f(14,i,j,nz+1); f(19,i,j,nz+1)=tmpval
         tmpval=f(16,i,j,nz+1); f(12,i,j,nz+1)=tmpval
         tmpval=f(18,i,j,nz+1); f(15,i,j,nz+1)=tmpval
         tmpval=f(20,i,j,nz+1); f(26,i,j,nz+1)=tmpval
         tmpval=f(23,i,j,nz+1); f(24,i,j,nz+1)=tmpval
         tmpval=f(25,i,j,nz+1); f(22,i,j,nz+1)=tmpval
         tmpval=f(27,i,j,nz+1); f(21,i,j,nz+1)=tmpval

         ix=i-cxs( 7); jy=j-cys( 7); tmpval=f( 7,ix,jy,nz); f( 7,i,jy+1,nz+1)=tmpval
         ix=i-cxs(13); jy=j-cys(13); tmpval=f(13,ix,jy,nz); f(13,i,jy+1,nz+1)=tmpval
         ix=i-cxs(14); jy=j-cys(14); tmpval=f(14,ix,jy,nz); f(14,i,jy+1,nz+1)=tmpval
         ix=i-cxs(16); jy=j-cys(16); tmpval=f(16,ix,jy,nz); f(16,i,jy+1,nz+1)=tmpval
         ix=i-cxs(18); jy=j-cys(18); tmpval=f(18,ix,jy,nz); f(18,i,jy+1,nz+1)=tmpval
         ix=i-cxs(20); jy=j-cys(20); tmpval=f(20,ix,jy,nz); f(20,i,jy+1,nz+1)=tmpval
         ix=i-cxs(23); jy=j-cys(23); tmpval=f(23,ix,jy,nz); f(23,i,jy+1,nz+1)=tmpval
         ix=i-cxs(25); jy=j-cys(25); tmpval=f(25,ix,jy,nz); f(25,i,jy+1,nz+1)=tmpval
         ix=i-cxs(27); jy=j-cys(27); tmpval=f(27,ix,jy,nz); f(27,i,jy+1,nz+1)=tmpval
      enddo
      enddo

!      do j = 0, ny+1
!      do i = 0, nx+1
!         f( 6,i,j,nz+1) =  f( 7,i-cxs( 7),j-cys( 7),nz)
!         f(17,i,j,nz+1) =  f(13,i-cxs(13),j-cys(13),nz)
!         f(19,i,j,nz+1) =  f(14,i-cxs(14),j-cys(14),nz)
!         f(12,i,j,nz+1) =  f(16,i-cxs(16),j-cys(16),nz)
!         f(15,i,j,nz+1) =  f(18,i-cxs(18),j-cys(18),nz)
!         f(26,i,j,nz+1) =  f(20,i-cxs(20),j-cys(20),nz)
!         f(24,i,j,nz+1) =  f(23,i-cxs(23),j-cys(23),nz)
!         f(22,i,j,nz+1) =  f(25,i-cxs(25),j-cys(25),nz)
!         f(21,i,j,nz+1) =  f(27,i-cxs(27),j-cys(27),nz)
!      enddo
!      enddo
   endif
!                           1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7
!   integer :: cxs(1:nl) = [0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1]
!   integer :: cys(1:nl) = [0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 0, 0,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1]
!   integer :: czs(1:nl) = [0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1]


! Periodic boundary conditions in i-direction.
   if (ibnd==0) then
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
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
!$cuf kernel do(2) <<<*,*>>>
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
!$cuf kernel do(2) <<<*,*>>>
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

   if ((jbnd < 10) .or. (kbnd < 10)) return


! === Edge along j=1 and k=1 === No slip
#ifdef CUDAS
!$cuf kernel do(2) <<<*,*>>>
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
!$cuf kernel do(2) <<<*,*>>>
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
!$cuf kernel do(2) <<<*,*>>>
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
!$cuf kernel do(2) <<<*,*>>>
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


   call cpufinish(icpu)

end subroutine
end module
