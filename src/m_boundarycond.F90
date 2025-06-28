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
   integer i,j,k,l,m,ja,ka
   integer, parameter :: icpu=8

   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions in i-direction
   if (ibnd==0) then
! Periodic boundary conditions in i-direction
      f(:,0   ,:,:)=f(:,nx,:,:)
      f(:,nx+1,:,:)=f(:,1 ,:,:)

   elseif (ibnd==1) then
! Inflow outflow boundary conditions in i-direction
      do k=1,nz
      do j=1,ny
         utmp(j,k)=uvel(k)
         vtmp(j,k)=0.0
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
            f(l,0,j,k)=f(l,0,j,k)-2.0*weights(l)*rtmp(ja,ka)*(cxs(l)*uvel(ka)+cys(l)*vtmp(ja,ka)+czs(l)*wtmp(ja,ka))/cs2
         enddo
         do l=2,nl-1,2
            tmp=f(l,0,j,k)
            if (cxs(l)==1)   f(l,0,j,k)=f(l+1,0,j,k)
            if (cxs(l+1)==1) f(l+1,0,j,k)=tmp
         enddo


!                           1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7
!   integer :: cxs(1:nl) = [0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1]
!   integer :: cys(1:nl) = [0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 0, 0,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1]
!   integer :: czs(1:nl) = [0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1]

      enddo
      enddo



      do k=0,nz+1
      do j=0,ny+1
         f(:,nx+1,j,k)=f(:,nx,j,k)
      enddo
      enddo

   elseif (ibnd==2) then
! Periodic boundary conditions with pressure gradient in i-direction
      f(:,0   ,:,:)=f(:,nx,:,:)
      f(:,nx+1,:,:)=f(:,1 ,:,:)
      call bndpressure(f,rho,u,v,w)
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Periodic boundary conditions in j-direction (sideways) periodic.
   if (jbnd==0) then
      f(:,:,0,:)   =f(:,:,ny,:)
      f(:,:,ny+1,:)=f(:,:,1,:)
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions in k-direction (vertical)
   if (kbnd==0) then ! Periodic boundary conditions in k-direction
      f(:,:,:,0)   =f(:,:,:,nz)
      f(:,:,:,nz+1)=f(:,: ,:,1)

   elseif (kbnd==1) then
! Outflow boundary on top
!      f(:,:,:,nz+1)=f(:,:,:,nz)
      do j = 1, ny
      do i = 1, nx
         f( 6,i,j,nz+1) = f( 7,i,j,nz)
         f(12,i,j,nz+1) = f(16,i,j,nz)
         f(15,i,j,nz+1) = f(18,i,j,nz)
         f(17,i,j,nz+1) = f(13,i,j,nz)
         f(19,i,j,nz+1) = f(14,i,j,nz)
         f(21,i,j,nz+1) = f(27,i,j,nz)
         f(22,i,j,nz+1) = f(25,i,j,nz)
         f(24,i,j,nz+1) = f(23,i,j,nz)
         f(26,i,j,nz+1) = f(20,i,j,nz)
!         do l = 1, nl
!            if (czs(l) < 0) then
!               ! Find specular reflection: (cx(j), cy(j), cz(j)) = (cx(i), cy(i), -cz(i))
!               do m = 1, nl
!                  if (cxs(m) == cxs(l) .and. cys(m) == cys(l) .and. czs(m) == -czs(l)) then
!                    f(l, i, j, nz+1) = f(m, i, j, nz)
!                    if ((j==ny/2).and.(i==nx/2)) print '(a,i2,a,i2,a)','f(',l,'i,j,nz+1) = f(',m,'i,j,nz)'
!                    exit
!                  endif
!               enddo
!            endif
!         enddo
      enddo
      enddo

   elseif (kbnd==2) then
! Outflow boundary on top and bottom
      f(:,:,:,0)   =f(:,:,:,1)
      f(:,:,:,nz+1)=f(:,:,:,nz)
   endif



   call cpufinish(icpu)

end subroutine
end module
