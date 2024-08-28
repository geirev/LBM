module m_actuatorline
contains
subroutine actuatorline(force,nx,ny,ipos,jpos,thetain,iradius,u,v,w)
   use m_readinfile, only : p2l,uini,turbrpm
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
   real theta                             ! work blade rotation angle
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

   real :: rho=1.0                        ! nondimensional density

   real, save :: CL(nrchords)             ! Lift coefficient for each chord along the blade
   real, save :: CD(nrchords)             ! Drag coefficient for each chord along the blade
   real :: forceL(nrchords)               ! Lift force along the blade
   real :: forceD(nrchords)               ! Drag force along the blade

   real :: gauss                          ! weigthing function for smeering force
   integer ichord                         ! chord counter
   integer, save :: ifirst=0
   real omeg

   real, save :: relma(nrchords),relmb(nrchords)
real vx,wx

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
      print *
      print '(a,17f8.2)','cl  :',cl(1:17)
      print '(a,17f8.2)','cd  :',cd(1:17)
      print '(a,17f8.2)','dc  :',dc(1:17)
      print '(a,17f8.2)','relm:',relm(1:17)
      print '(a,17f8.2)','ra  :',relma(1:17)
      print '(a,17f8.2)','chor:',chord(1:17)
   endif

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


! Rotation speed in radians/s
   omega=pi2*turbrpm/60.0

   if (ifirst==1) then
      print *,'Omega  (RPM)                 =',turbrpm
      print *,'Omega  (radians/s)           =',omega
! Tipspeed(Tip speed can be determined from the rotational speed, which is ωR where ω is the rotational
! speed in radians per second and R is the radius of the turbine in meters.)
      print *,'tipspeed R*Omega m/s         =',real(iradius)*p2l%length*omega
! Tipspeed ratio 8 m/s winds (For three blades, a TSR of 6 to 7 is optimal. If it is less, not all
! the available energy is captured; if it is more, the blades move into an area of turbulence from the last
! blade and are not as efficient.
      print *,'tipspeed ratio R*Omega/uini  =',real(iradius)*p2l%length*omega/(uini*p2l%vel)
   endif

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

! Here we could replace this with bilinear or cubic interpolation, to find u(xb,yb), v(xb,yb), w(xb,yb)
! in the square  ic=int(xb), ic+1,  jc=int(yb), jc+1
! Non-dimensional velocity components from model
         ux = u(ic,jc)
         vx = v(ic,jc)
         wx = w(ic,jc)

! Non-dimensional utheta
         utheta =  omega*sqrt((xb-x0)**2+(yb-y0)**2+0.000001) - wx*cos(theta) - vx*sin(theta)

! Non-dimensional urel**2
         urel2  = ux**2 + utheta**2

! Non-dimensional dynamic pressure
         q= 0.5 * rho * urel2

! Non-dimensional lift and drag forces per blade element
         forceL(ichord) = q * chord(ichord) * dc(ichord) * CL(ichord)
         forceD(ichord) = q * chord(ichord) * dc(ichord) * CD(ichord)

! Non-dimensional lift and drag forces per blade element unit length
         forceL(ichord) = forceL(ichord)/dc(ichord)
         forceD(ichord) = forceD(ichord)/dc(ichord)

! Flowangle between relative windspeed and rotor plane
         phi    = atan2(ux,utheta)    ! pi/2.0 - atan2(ux,utheta)

         if (ifirst==1) then
            if (ichord == 1) then
               print *,'Blade :',iblade
               print '(2a8,2a3,9a8)','  ichord','    dist',' ic',' jc','      ux','      vx','      wx','  utheta',&
                                     '     phi','    urel','      dc','  forceL','  forceD'
            endif
            print '(i8.8,f8.4,2I3,9f8.4)',ichord,sqrt((xb-x0)**2+(yb-y0)**2+0.0001),ic,jc,ux,vx,wx,utheta,&
                                          phi*360.0/pi2,sqrt(urel2),dc(ichord),forceL(ichord),forceD(ichord)

         endif
      enddo



! Computing the forces for gridpoints located at (x,y)
      do j=ja,jb
         y=real(j)
         do i=ia,ib
            x=real(i)

            t=min(1.0, ((x-x0)*a + (y-y0)*b)/radius**2)
            t=((x-x0)*a + (y-y0)*b)/radius**2

            ! location on blade
            xp=x0+t*a
            yp=y0+t*b

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

! Note the negative v, is because of coordinate system orientation v to left and w up in rotor plan
            if ((0.01 < t).and.(t <= 1.0) .and. (nc > 0)) then
               gauss=(1.0/(sqrt(pi)*eps)) * exp(-((x-xp)**2 + (y-yp)**2)/eps**2)
               force(i,j,1)=force(i,j,1)+(forceL(nc)*cos(phi) + forceD(nc)*sin(phi))*gauss
               force(i,j,2)=force(i,j,2)-(forceL(nc)*sin(phi) - forceD(nc)*cos(phi))*sin(theta)*gauss
               force(i,j,3)=force(i,j,3)+(forceL(nc)*sin(phi) - forceD(nc)*cos(phi))*cos(theta)*gauss
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
