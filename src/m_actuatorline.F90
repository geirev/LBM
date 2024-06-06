module m_actuatorline
contains
subroutine actuatorline(force,nx,ny,ipos,jpos,thetain,omegain,u,v,w)
   use m_readinfile, only : p2l
   use mod_nrl5mw
   implicit none
   integer, intent(in)    :: nx           ! dimension of horizontal axis of turbine plane
   integer, intent(in)    :: ny           ! dimension of vertical   axis of turbine plane
   integer, intent(in)    :: ipos         ! horizontal gridpoint for center of turbine
   integer, intent(in)    :: jpos         ! vertical   gridpoint for center of turbine


   real,    intent(in)    :: u(nx,ny)     ! u velocity through the turbine section
   real,    intent(in)    :: v(nx,ny)     ! v velocity at the turbine section
   real,    intent(in)    :: w(nx,ny)     ! w velocity at the turbine section
   real,    intent(in)    :: thetain      ! input rotation angle of blade one
   real,    intent(in)    :: omegain      ! rotation speed in RPM
   real,    intent(inout) :: force(nx,ny,3) ! actuator line force added to force



   real, parameter :: dx=1.0
   real, parameter :: pi=3.1415926535
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0
   real, parameter :: eps=1.6             ! smoothing distance in Gaussian kernel

   real x,y                               ! x-y locatiuon along a blade
   real xy                                ! length along the blade from x1,y1 to xp,yp
   real x0,y0                             ! x-y locatiuon of a grid point
   real x1,y1                             ! x-y location of turbine center
   real x2,y2                             ! x-y location of blade tip
   real xp,yp                             ! x-y locatiopn of point on blade closest to x0-y0
   real a,b,t                             ! work variables
   real theta                             ! work rotation angle
   real omega                             ! Rotation speed in radians per second
   real phi                               ! Vrel angle
   real q                                 ! Dynamic pressure
   real ux                                ! inflow velocity
   real utheta                            ! tangential velocity from v and w
   real urel2                             ! relative velocity squared
   real, save :: radius                            ! turbine radius in grid cells
   integer, save :: iradius                        ! number of gridpoints for blade length
   integer i,ia,ib                        ! counter and loop limits in horizontal direction
   integer j,ja,jb                        ! counter and loop limits in vertical direction
   integer ic,jc,nc                       ! gridpoint along a blade
   integer iblade                         ! blade counter

   real :: rho=1.0                        ! nondimentional density

   real, save :: CL(nrchords)                   ! Lift coefficient for each chord along the blade
   real, save :: CD(nrchords)                   ! Drag coefficient for each chord along the blade
   real :: forceL(nrchords)               ! Lift force along the blade
   real :: forceD(nrchords)               ! Drag force along the blade

   real :: gauss                          ! weigthing function for smeering force
   integer ichord                         ! chord counter
   integer, save :: ifirst=0

   real, save :: relma(nrchords),relmb(nrchords)

   ifirst=ifirst+1
   force=0.0


   if (ifirst == 1) then
! Initialize lift and drag coefficients
      CL(:)=0.2
      CD(:)=0.2

! Nondimensional turbine parameters
      relm=relm/p2l%length
      dc=dc/p2l%length
      chord=chord/p2l%length
      rotorradius=rotorradius/p2l%length
      hubradius=hubradius/p2l%length

      do ichord=1,nrchords
         relma(ichord)=relm(ichord)-dc(ichord)*0.5
         relmb(ichord)=relm(ichord)+dc(ichord)*0.5
      enddo

! Blade length in number of gridpoints
      radius=rotorradius+hubradius     ! radius comes from mod_nrl5Mw
      iradius=nint(radius)             ! radius in number of gridcells
   endif
!   print '(a,10f8.2)','dc  :',dc(1:10)
!   print '(a,10f8.2)','relm:',relm(1:10)
!   print '(a,10f8.2)','ra  :',relma(1:10)
!   print '(a,10f8.2)','chor:',chord(1:10)

! Set rotation angle
   theta=thetain

! Center point of turbine
   x1=real(ipos)
   y1=real(jpos)

! Loop limits for computing force
   ia=ipos-iradius-2
   ib=ipos+iradius+2
   ja=jpos-iradius-2
   jb=jpos+iradius+2

   if (ifirst==1) then
      open(21,file='dist.dat',status='unknown')
      write(21,*)'TITLE = "Dist"'
      write(21,*)'VARIABLES = "i-index" "j-index" "u" "v" "w" "u2"'
      close(21)
   endif


! rotation speed in radians/s
   omega=pi2*omegain/60.0

! nondimesonal omega
   omega=omega*p2l%time

! Running through blades
   do iblade=1,3

! Find parameterization of line from center (x1,y1) to blade tip (x2,y2)
      x2=x1+radius*cos(theta)
      y2=y1+radius*sin(theta)

      a=x2-x1
      b=y2-y1

! Compute lift and drag force along blade
      do ichord=1,nrchords
         ! Finding gridpoint index closest to chord number ichord to extract velocity
         x=x1+relm(ichord)*cos(theta)
         y=y1+relm(ichord)*sin(theta)
         ic=nint(x)
         jc=nint(y)

         ux     = u(ic,jc)
         utheta = sqrt( omega*sqrt((x-x1)**2+(y-y1)**2+0.001) - v(ic,jc)*cos(theta) + w(ic,jc)*sin(theta) )
         phi    = 0.5*pi-atan2(ux,utheta)
         urel2  = ux**2 + utheta**2

! Dynamic pressure
         q= 0.5 * rho * urel2

! Lift and drag forces per element
         forceL(ichord) = q * chord(ichord) * dc(ichord) * CL(ichord)
         forceD(ichord) = q * chord(ichord) * dc(ichord) * CD(ichord)

         forceL(ichord) = forceL(ichord)/dc(ichord)
         forceD(ichord) = forceD(ichord)/dc(ichord)

!         print '(3(a,i3,tr1),7(a,f10.3,tr1))','ich=',ichord,'ic=',ic,'jc=',jc,'ux=',ux,'ut=',utheta,&
!                                'theta=',theta*360.0/pi2,'phi=',phi*360.0/pi2,'fL=',forceL(ichord),'fF=',forceD(ichord)
      enddo
!      print '(a,10f8.2)','forcL:',forceL(1:10)
!      print '(a,10f8.2)','forcD:',forceD(1:10)


!      if (ifirst < 5) then
!         open(10,file='nrl.dat',position='append')
!            do i=1,nrchords
!               write(10,'(i3,6f12.2)')i,relm(i),hubradius+sum(dc(1:i-1))+dc(i)/2.0,dc(i),chord(i),forceL(i),forceD(i)
!            enddo
!         close(10)
!      endif

! Computing the forces for gripoints located at (x0,y0)
      do j=ja,jb
         y0=real(j)
         do i=ia,ib
            x0=real(i)

            t=min(1.0, ((x0-x1)*a + (y0-y1)*b)/radius**2)

            ! location on blade
            xp=x1+t*a
            yp=y1+t*b

            !nc=nint(t*nrchords)
            xy=sqrt((xp-x1)**2 + (yp-y1)**2)    ! length from (x1,y1) to location (xp,yp) on blade
            nc=0
            do ichord=1,nrchords
               if (relma(ichord) <= xy .and. xy < relmb(ichord)) then
                  nc=ichord
                  exit
               endif
            enddo
            if (xy >= hubradius .and. nc==0) then
                  nc=nrchords
            endif

!            if (nc == 0) print *,'problem (nc): ',nc,t,i,j,xy,hubradius


            if ((0.01 < t).and.(t <= 1.0) .and. (nc > 0)) then
               gauss=(1.0/(sqrt(pi)*eps)) * exp(-((xp-x0)**2 + (yp-y0)**2)/eps**2)
               force(i,j,1)=force(i,j,1)+(forceL(nc)*cos(phi) + forceD(nc)*sin(phi))*gauss
               force(i,j,3)=force(i,j,3)+(forceL(nc)*sin(phi) - forceD(nc)*cos(phi))*cos(theta)*gauss
               force(i,j,2)=force(i,j,2)-(forceL(nc)*sin(phi) - forceD(nc)*cos(phi))*sin(theta)*gauss
            endif

         enddo
      enddo

      theta=theta+rad120
   enddo

   if (ifirst <= 200) then
      open(21,file='dist.dat',position='append')
      write(21,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',nx,', J=',ny,', K=1'
      write(21,'(30I5)')((i,i=1,nx),j=1,ny)
      write(21,'(30I5)')((j,i=1,nx),j=1,ny)
      write(21,'(10(1x,e12.5))')((force(i,j,1),i=1,nx),j=1,ny)
      write(21,'(10(1x,e12.5))')((force(i,j,2),i=1,nx),j=1,ny)
      write(21,'(10(1x,e12.5))')((force(i,j,3),i=1,nx),j=1,ny)
      write(21,'(10(1x,e12.5))')((sqrt(force(i,j,1)**2+force(i,j,2)**2+force(i,j,3)**2),i=1,nx),j=1,ny)
      close(21)
   endif

end subroutine
end module
