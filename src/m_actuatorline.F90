module m_actuatorline
contains
subroutine actuatorline(force,nx,ny,ipos,jpos,thetain,omegain,iradius,u,v,w)
   use m_readinfile, only : p2l
   use mod_nrl5mw
   use m_readfoildata
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
   integer, intent(inout) :: iradius      ! number of gridpoints for blade length computed at first call
   real,    intent(inout) :: force(nx,ny,3) ! actuator line force added to force



   real, parameter :: dx=1.0
   real, parameter :: pi=3.1415926535
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0
   real, parameter :: eps=1.6             ! smoothing distance in Gaussian kernel

   real x0,y0                             ! x-y location of turbine center
   real xb,yb                             ! x-y locatiuon along a blade
   real x2,y2                             ! x-y location of blade tip
   real x,y                               ! x-y locatiuon of a grid point
   real xp,yp                             ! x-y locatiopn of point on blade closest to x-y
   real xy                                ! length along the blade from x0,y0 to xp,yp
   real a,b,t                             ! work variables
   real theta                             ! work rotation angle
   real omega                             ! Rotation speed in radians per second
   real phi                               ! Vrel angle
   real q                                 ! Dynamic pressure
   real ux                                ! inflow velocity
   real utheta                            ! tangential velocity from v and w
   real urel2                             ! relative velocity squared
   real, save :: radius                   ! turbine radius in grid cells
   integer i,ia,ib                        ! counter and loop limits in horizontal direction
   integer j,ja,jb                        ! counter and loop limits in vertical direction
   integer ic,jc,nc                       ! gridpoint along a blade
   integer iblade                         ! blade counter

   real :: rho=1.0                        ! nondimentional density

   real, save :: CL(nrchords)             ! Lift coefficient for each chord along the blade
   real, save :: CD(nrchords)             ! Drag coefficient for each chord along the blade
   real :: forceL(nrchords)               ! Lift force along the blade
   real :: forceD(nrchords)               ! Drag force along the blade

   real :: gauss                          ! weigthing function for smeering force
   integer ichord                         ! chord counter
   integer, save :: ifirst=0
   real omeg

   real, save :: relma(nrchords),relmb(nrchords)

   ifirst=ifirst+1
   force=0.0


   if (ifirst == 1) then
! Read lift and drag from foil data files
      call readfoildata(cl,cd)

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
      print *,'iradius=',iradius
      print *,'iD=',2*iradius
      print *,'3*iD  =',3*2*iradius
      print *,'29*iD =',29*2*iradius
      print *,'omega(12.1 RPM/60 s)=',12.1/60.0
      print *,'tipspeed R*Omega    =',63.0*12.1/60.0
      print *,'tipspeed ratio      =',(63.0*12.1/60.0)/8.0
      print *
   endif
!   print '(a,10f8.2)','dc  :',dc(1:10)
!   print '(a,10f8.2)','relm:',relm(1:10)
!   print '(a,10f8.2)','ra  :',relma(1:10)
!   print '(a,10f8.2)','chor:',chord(1:10)

! Set rotation angle
   theta=thetain

! Center point of turbine
   x0=real(ipos)
   y0=real(jpos)

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

! Find parameterization of line from center (x0,y0) to blade tip (x2,y2)
      x2=x0+radius*cos(theta)
      y2=y0+radius*sin(theta)

      a=x2-x0
      b=y2-y0

! Compute lift and drag force along blade
      do ichord=1,nrchords
         ! Finding gridpoint index closest to chord number ichord to extract velocity
         xb=x0+relm(ichord)*cos(theta)
         yb=y0+relm(ichord)*sin(theta)
         ic=nint(xb)
         jc=nint(yb)

         ux     = u(ic,jc)
         utheta = sqrt( omega*sqrt((xb-x0)**2+(yb-y0)**2+0.001) + v(ic,jc)*cos(theta) + w(ic,jc)*sin(theta) )
         phi    = 0.5*pi-atan2(ux,utheta)
         urel2  = ux**2 + utheta**2

! Dynamic pressure
         q= 0.5 * rho * urel2

! Lift and drag forces per element
         forceL(ichord) = q * chord(ichord) * dc(ichord) * CL(ichord)
         forceD(ichord) = q * chord(ichord) * dc(ichord) * CD(ichord)

         forceL(ichord) = forceL(ichord)/dc(ichord)
         forceD(ichord) = forceD(ichord)/dc(ichord)
      enddo


! Computing the forces for gripoints located at (x,y)
      do j=ja,jb
         y=real(j)
         do i=ia,ib
            x=real(i)

            t=min(1.0, ((x-x0)*a + (y-y0)*b)/radius**2)

            ! location on blade
            xp=x0+t*a
            yp=y0+t*b

            !nc=nint(t*nrchords)
            xy=sqrt((xp-x0)**2 + (yp-y0)**2)    ! length from (x0,y0) to location (xp,yp) on blade
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

            if ((0.01 < t).and.(t <= 1.0) .and. (nc > 0)) then
               gauss=(1.0/(sqrt(pi)*eps)) * exp(-((xp-x)**2 + (yp-y)**2)/eps**2)
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
