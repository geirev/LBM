module m_boundarycond
contains
subroutine boundarycond(f,rho,u,v,w,uvel)
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_bndpressure
   use m_fequilscalar
   use m_wtime
   implicit none
   real, intent(inout):: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)   :: rho(nx,ny,nz)
   real, intent(in)   :: u(nx,ny,nz)
   real, intent(in)   :: v(nx,ny,nz)
   real, intent(in)   :: w(nx,ny,nz)
   real, intent(in)   :: uvel(nz)
   real :: rtmp(ny,nz)
   real :: utmp(ny,nz)
   real :: vtmp(ny,nz)
   real :: wtmp(ny,nz)
!   real uin,rhoin
   real tmp
   integer i,j,k,l,m,ja,ka,ip
   integer, parameter :: icpu=8
!   integer bouncefreey(nl)
!   integer bouncefreez(nl)

   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Closed boundary conditions in i-direction

   if (ibnd==1) then
! Inflow outflow boundary conditions in i-direction
      do k=1,nz
      do j=1,ny
         utmp(j,k)=uvel(k)
         vtmp(j,k)=0.5*uvel(k)
         wtmp(j,k)=0.0
         rtmp(j,k)=rho0 !rho(1,j,k)
      enddo
      enddo

      do k=0,nz+1
      ka=min(max(k,1),nz)
      do j=0,ny+1
         ja=min(max(j,1),ny)
!         f(1:nl,0,j,k)=fequilscalar(rtmp(ja,ka),utmp(ja,ka),vtmp(ja,ka),wtmp(ja,ka))
         f(1:nl,0,j,k)=f(1:nl,1,j,k)
!         rhoin=0.0
!         do l=1,nl
!            if (cxs(l) <= 0) rhoin= rhoin + f(l,0,ja,ka)
!         enddo
!         uin=utmp(ja,ka)
!         rhoin = rhoin / (1.0 - uin)

         do l=1,nl
            f(l,0,j,k)=f(l,0,j,k)-2.0*weights(l)*rtmp(ja,ka)*(cxs(l)*utmp(ja,ka)+cys(l)*vtmp(ja,ka)+czs(l)*wtmp(ja,ka))/cs2
         enddo
         do l=2,nl-1,2
            tmp=f(l,0,j,k)
            if (cxs(l)==1)   f(l,0,j,k)=f(l+1,0,j,k)
            if (cxs(l+1)==1) f(l+1,0,j,k)=tmp
         enddo



      enddo
      enddo



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
      do k=1,nz
      do i=1,nx
         f( 4,i,0,k) = f( 5,i,0,k)
         f( 8,i,0,k) = f( 9,i,0,k)
         f(11,i,0,k) = f(10,i,0,k)
         f(14,i,0,k) = f(15,i,0,k)
         f(19,i,0,k) = f(18,i,0,k)
         f(20,i,0,k) = f(21,i,0,k)
         f(23,i,0,k) = f(22,i,0,k)
         f(24,i,0,k) = f(25,i,0,k)
         f(26,i,0,k) = f(27,i,0,k)

         f( 5,i,0,k) = f( 5,i-cxs( 5),1,k-czs( 5))
         f( 9,i,0,k) = f( 9,i-cxs( 9),1,k-czs( 9))
         f(10,i,0,k) = f(10,i-cxs(10),1,k-czs(10))
         f(15,i,0,k) = f(15,i-cxs(15),1,k-czs(15))
         f(18,i,0,k) = f(18,i-cxs(18),1,k-czs(18))
         f(21,i,0,k) = f(21,i-cxs(21),1,k-czs(21))
         f(22,i,0,k) = f(22,i-cxs(22),1,k-czs(22))
         f(25,i,0,k) = f(25,i-cxs(25),1,k-czs(25))
         f(27,i,0,k) = f(27,i-cxs(27),1,k-czs(27))
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
      do k=1,nz
      do i=1,nx
         f( 4,i,0,k) = f( 5,i,0,k)
         f(11,i,0,k) = f( 9,i,0,k)
         f( 8,i,0,k) = f(10,i,0,k)
         f(19,i,0,k) = f(15,i,0,k)
         f(14,i,0,k) = f(18,i,0,k)
         f(24,i,0,k) = f(21,i,0,k)
         f(26,i,0,k) = f(22,i,0,k)
         f(20,i,0,k) = f(25,i,0,k)
         f(23,i,0,k) = f(27,i,0,k)

         f( 5,i,0,k) = f( 5,i-cxs( 5),1,k-czs( 5))
         f( 9,i,0,k) = f( 9,i-cxs( 9),1,k-czs( 9))
         f(10,i,0,k) = f(10,i-cxs(10),1,k-czs(10))
         f(15,i,0,k) = f(15,i-cxs(15),1,k-czs(15))
         f(18,i,0,k) = f(18,i-cxs(18),1,k-czs(18))
         f(21,i,0,k) = f(21,i-cxs(21),1,k-czs(21))
         f(22,i,0,k) = f(22,i-cxs(22),1,k-czs(22))
         f(25,i,0,k) = f(25,i-cxs(25),1,k-czs(25))
         f(27,i,0,k) = f(27,i-cxs(27),1,k-czs(27))
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
      do k=1,nz
      do i=1,nx
         f( 5,i,ny+1,k) = f( 4,i,ny+1,k)
         f( 9,i,ny+1,k) = f( 8,i,ny+1,k)
         f(10,i,ny+1,k) = f(11,i,ny+1,k)
         f(15,i,ny+1,k) = f(14,i,ny+1,k)
         f(18,i,ny+1,k) = f(19,i,ny+1,k)
         f(21,i,ny+1,k) = f(20,i,ny+1,k)
         f(22,i,ny+1,k) = f(23,i,ny+1,k)
         f(25,i,ny+1,k) = f(24,i,ny+1,k)
         f(27,i,ny+1,k) = f(26,i,ny+1,k)

         f( 4,i,ny+1,k) = f( 4,i-cxs( 4),ny,k-czs( 4))
         f( 8,i,ny+1,k) = f( 8,i-cxs( 8),ny,k-czs( 8))
         f(11,i,ny+1,k) = f(11,i-cxs(11),ny,k-czs(11))
         f(14,i,ny+1,k) = f(14,i-cxs(14),ny,k-czs(14))
         f(19,i,ny+1,k) = f(19,i-cxs(19),ny,k-czs(19))
         f(20,i,ny+1,k) = f(20,i-cxs(20),ny,k-czs(20))
         f(23,i,ny+1,k) = f(23,i-cxs(23),ny,k-czs(23))
         f(24,i,ny+1,k) = f(24,i-cxs(24),ny,k-czs(24))
         f(26,i,ny+1,k) = f(26,i-cxs(26),ny,k-czs(26))
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
      do k=1,nz
      do i=1,nx
         f( 5,i,ny+1,k) = f( 4,i,ny+1,k)
         f(10,i,ny+1,k) = f( 8,i,ny+1,k)
         f( 9,i,ny+1,k) = f(11,i,ny+1,k)
         f(18,i,ny+1,k) = f(14,i,ny+1,k)
         f(15,i,ny+1,k) = f(19,i,ny+1,k)
         f(25,i,ny+1,k) = f(20,i,ny+1,k)
         f(27,i,ny+1,k) = f(23,i,ny+1,k)
         f(21,i,ny+1,k) = f(24,i,ny+1,k)
         f(22,i,ny+1,k) = f(26,i,ny+1,k)

         f( 4,i,ny+1,k) = f( 4,i-cxs( 4),ny,k-czs( 4))
         f( 8,i,ny+1,k) = f( 8,i-cxs( 8),ny,k-czs( 8))
         f(11,i,ny+1,k) = f(11,i-cxs(11),ny,k-czs(11))
         f(14,i,ny+1,k) = f(14,i-cxs(14),ny,k-czs(14))
         f(19,i,ny+1,k) = f(19,i-cxs(19),ny,k-czs(19))
         f(20,i,ny+1,k) = f(20,i-cxs(20),ny,k-czs(20))
         f(23,i,ny+1,k) = f(23,i-cxs(23),ny,k-czs(23))
         f(24,i,ny+1,k) = f(24,i-cxs(24),ny,k-czs(24))
         f(26,i,ny+1,k) = f(26,i-cxs(26),ny,k-czs(26))
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
      do j=1,ny
      do i=1,nx
         f( 7,i,j,0) = f( 6,i,j,0)
         f(13,i,j,0) = f(12,i,j,0)
         f(14,i,j,0) = f(15,i,j,0)
         f(16,i,j,0) = f(17,i,j,0)
         f(18,i,j,0) = f(19,i,j,0)
         f(20,i,j,0) = f(21,i,j,0)
         f(23,i,j,0) = f(22,i,j,0)
         f(25,i,j,0) = f(24,i,j,0)
         f(27,i,j,0) = f(26,i,j,0)

         f( 6,i,j,0) = f( 6,i-cxs( 6),j-cys( 6),1)
         f(12,i,j,0) = f(12,i-cxs(12),j-cys(12),1)
         f(15,i,j,0) = f(15,i-cxs(15),j-cys(15),1)
         f(17,i,j,0) = f(17,i-cxs(17),j-cys(17),1)
         f(19,i,j,0) = f(19,i-cxs(19),j-cys(19),1)
         f(21,i,j,0) = f(21,i-cxs(21),j-cys(21),1)
         f(22,i,j,0) = f(22,i-cxs(22),j-cys(22),1)
         f(24,i,j,0) = f(24,i-cxs(24),j-cys(24),1)
         f(26,i,j,0) = f(26,i-cxs(26),j-cys(26),1)
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
      do j=1,ny
      do i=1,nx
         f( 7,i,j,0) = f( 6,i,j,0)
         f(16,i,j,0) = f(12,i,j,0)
         f(18,i,j,0) = f(15,i,j,0)
         f(13,i,j,0) = f(17,i,j,0)
         f(14,i,j,0) = f(19,i,j,0)
         f(27,i,j,0) = f(21,i,j,0)
         f(25,i,j,0) = f(22,i,j,0)
         f(23,i,j,0) = f(24,i,j,0)
         f(20,i,j,0) = f(26,i,j,0)

         f( 6,i,j,0) = f( 6,i-cxs( 6),j-cys( 6),1)
         f(12,i,j,0) = f(12,i-cxs(12),j-cys(12),1)
         f(15,i,j,0) = f(15,i-cxs(15),j-cys(15),1)
         f(17,i,j,0) = f(17,i-cxs(17),j-cys(17),1)
         f(19,i,j,0) = f(19,i-cxs(19),j-cys(19),1)
         f(21,i,j,0) = f(21,i-cxs(21),j-cys(21),1)
         f(22,i,j,0) = f(22,i-cxs(22),j-cys(22),1)
         f(24,i,j,0) = f(24,i-cxs(24),j-cys(24),1)
         f(26,i,j,0) = f(26,i-cxs(26),j-cys(26),1)
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
      do j=1,ny+1
      do i=1,nx+1
         f( 6,i,j,nz+1) = f( 7,i,j,nz+1)
         f(12,i,j,nz+1) = f(13,i,j,nz+1)
         f(15,i,j,nz+1) = f(14,i,j,nz+1)
         f(17,i,j,nz+1) = f(16,i,j,nz+1)
         f(19,i,j,nz+1) = f(18,i,j,nz+1)
         f(21,i,j,nz+1) = f(20,i,j,nz+1)
         f(22,i,j,nz+1) = f(23,i,j,nz+1)
         f(24,i,j,nz+1) = f(25,i,j,nz+1)
         f(26,i,j,nz+1) = f(27,i,j,nz+1)

         f( 7,i,j,nz+1) = f( 7,i-cxs( 7),j-cys( 7),nz)
         f(13,i,j,nz+1) = f(13,i-cxs(13),j-cys(13),nz)
         f(14,i,j,nz+1) = f(14,i-cxs(14),j-cys(14),nz)
         f(16,i,j,nz+1) = f(16,i-cxs(16),j-cys(16),nz)
         f(18,i,j,nz+1) = f(18,i-cxs(18),j-cys(18),nz)
         f(20,i,j,nz+1) = f(20,i-cxs(20),j-cys(20),nz)
         f(23,i,j,nz+1) = f(23,i-cxs(23),j-cys(23),nz)
         f(25,i,j,nz+1) = f(25,i-cxs(25),j-cys(25),nz)
         f(27,i,j,nz+1) = f(27,i-cxs(27),j-cys(27),nz)
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

      do j = 1, ny
      do i = 1, nx
         f( 6,i,j,nz+1) =  f( 7,i,j,nz+1)
         f(17,i,j,nz+1) =  f(13,i,j,nz+1)
         f(19,i,j,nz+1) =  f(14,i,j,nz+1)
         f(12,i,j,nz+1) =  f(16,i,j,nz+1)
         f(15,i,j,nz+1) =  f(18,i,j,nz+1)
         f(26,i,j,nz+1) =  f(20,i,j,nz+1)
         f(24,i,j,nz+1) =  f(23,i,j,nz+1)
         f(22,i,j,nz+1) =  f(25,i,j,nz+1)
         f(21,i,j,nz+1) =  f(27,i,j,nz+1)

         f( 7,i,j,nz+1) =  f( 7,i-cxs( 7),j-cys( 7),nz)
         f(13,i,j,nz+1) =  f(13,i-cxs(13),j-cys(13),nz)
         f(14,i,j,nz+1) =  f(14,i-cxs(14),j-cys(14),nz)
         f(16,i,j,nz+1) =  f(16,i-cxs(16),j-cys(16),nz)
         f(18,i,j,nz+1) =  f(18,i-cxs(18),j-cys(18),nz)
         f(20,i,j,nz+1) =  f(20,i-cxs(20),j-cys(20),nz)
         f(23,i,j,nz+1) =  f(23,i-cxs(23),j-cys(23),nz)
         f(25,i,j,nz+1) =  f(25,i-cxs(25),j-cys(25),nz)
         f(27,i,j,nz+1) =  f(27,i-cxs(27),j-cys(27),nz)
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
      f(:,0   ,:,:)=f(:,nx,:,:)
      f(:,nx+1,:,:)=f(:,1 ,:,:)
   endif

! Periodic boundary conditions in j-direction.
   if (jbnd==0) then
      f(:,:,0,:)   =f(:,:,ny,:)
      f(:,:,ny+1,:)=f(:,:,1,:)
   endif

! Periodic boundary conditions in k-direction.
   if (kbnd==0) then
      f(:,:,:,0)   =f(:,:,:,nz)
      f(:,:,:,nz+1)=f(:,: ,:,1)
   endif





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bounce-back edges in case of closed boundaries in both y and z directions
! Only needed in case of closed boundaries in both y and z directions

! The free-slip edges have not yet been implemented and I didn't care about the corners yet.

   if ((jbnd < 10) .or. (kbnd < 10)) return
   k=1; j=1; i=1; l=1; m=1; ip=1


! === Edge along j=1 and k=1 === No slip
    do i=1,nx
       do l=1,nl
          ip = bounce(l)
          if (cys(l) < 0.0 .and. czs(l) < 0.0) then
            f(ip,i,0,0) = f(l,i,0,0)
            f(l, i,0,0) = f(l, i - cxs(l), 0 - cys(l), 0 - czs(l))
          endif
          if (cys(l) <= 0.0 .and. czs(l) < 0.0) then
             f(ip,i,1,0) = f(l,i,1,0)
             f( l,i,1,0) = f(l, i - cxs(l), 1 - cys(l), 0 - czs(l))
          endif
          if (cys(l) < 0.0 .and. czs(l) <= 0.0) then
             f(ip,i,0,1) = f(l,i,0,1)
             f( l,i,0,1) = f(l,i - cxs(l), 0 - cys(l), 1 - czs(l))
          endif
       enddo
    enddo

! === Edge along j=ny and k=1 === No slip
    do i=1,nx
       do l=1,nl
          ip = bounce(l)
          if (cys(l) > 0.0 .and. czs(l) < 0.0) then
            f(ip,i,ny+1,0) = f(l,i,ny+1,0)
            f(l, i,ny+1,0) = f(l, i - cxs(l), ny+1 - cys(l), 0 - czs(l))
          endif
          if (cys(l) >= 0.0 .and. czs(l) > 0.0) then
             f(ip,i,ny,0) = f(l,i,ny,0)
             f( l,i,ny,0) = f(l, i - cxs(l), ny - cys(l), 0 - czs(l))
          endif
          if (cys(l) > 0.0 .and. czs(l) >= 0.0) then
             f(ip,i,ny+1,1) = f(l,i,ny+1,1)
             f( l,i,ny+1,1) = f(l,i - cxs(l), ny+1 - cys(l), 1 - czs(l))
          endif
       enddo
    enddo

! === Edge along j=1 and k=nz === No slip
    do i=1,nx
       do l=1,nl
          ip = bounce(l)
          if (cys(l) < 0.0 .and. czs(l) > 0.0) then
            f(ip,i,0,nz+1) = f(l,i,0,nz+1)
            f(l, i,0,nz+1) = f(l, i - cxs(l), 0 - cys(l), nz+1 - czs(l))
          endif
          if (cys(l) <= 0.0 .and. czs(l) > 0.0) then
             f(ip,i,1,nz+1) = f(l,i,1,nz+1)
             f( l,i,1,nz+1) = f(l,i - cxs(l), 1 - cys(l), nz+1 - czs(l))
          endif
          if (cys(l) < 0.0 .and. czs(l) >= 0.0) then
             f(ip,i,0,nz) = f(l,i,0,nz)
             f( l,i,0,nz) = f(l,i - cxs(l), 0 - cys(l), nz - czs(l))
          endif
       enddo
    enddo

! === Edge along j=ny and k=nz === No slip
    do i=1,nx
       do l=1,nl
          ip = bounce(l)
          if (cys(l) > 0.0 .and. czs(l) > 0.0) then
            f(ip,i,ny+1,nz+1) = f(l,i,ny+1,nz+1)
            f(l, i,ny+1,nz+1) = f(l, i - cxs(l), ny+1 - cys(l), nz+1 - czs(l))
          endif
          if (cys(l) >= 0.0 .and. czs(l) > 0.0) then
             f(ip,i,ny,nz+1) = f(l,i,ny,nz+1)
             f( l,i,ny,nz+1) = f(l, i - cxs(l), ny - cys(l), nz+1 - czs(l))
          endif
          if (cys(l) > 0.0 .and. czs(l) >= 0.0) then
             f(ip,i,ny+1,nz) = f(l,i,ny+1,nz)
             f( l,i,ny+1,nz) = f(l,i - cxs(l), ny+1 - cys(l), nz - czs(l))
          endif
       enddo
    enddo


   call cpufinish(icpu)

end subroutine
end module
