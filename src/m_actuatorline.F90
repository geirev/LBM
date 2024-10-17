module m_actuatorline
contains
subroutine actuatorline(force,nx,ny,ipos,jpos,thetain,iradius,u,v,w,rho,ieps)
   use m_readinfile, only : p2l,uini,turbrpm,nturbines,pitchangle,itiploss,tipspeedratio
   use mod_nrel5mw
   use m_nrelreadfoil
   use m_nrelliftdrag
   use m_bilin
   implicit none
   integer, intent(in)    :: nx           ! dimension of horizontal axis of turbine plane
   integer, intent(in)    :: ny           ! dimension of vertical   axis of turbine plane
   integer, intent(in)    :: ipos         ! horizontal gridpoint for center of turbine
   integer, intent(in)    :: jpos         ! vertical   gridpoint for center of turbine
   integer, intent(in)    :: ieps         ! Number of gridpoints to extend forcing orthogonal to rotor plane


   real,    intent(in)    :: u(nx,ny)     ! u velocity through the turbine section
   real,    intent(in)    :: v(nx,ny)     ! v velocity at the turbine section
   real,    intent(in)    :: w(nx,ny)     ! w velocity at the turbine section
   real,    intent(in)    :: rho(nx,ny)   ! density at the turbine section
   real,    intent(in)    :: thetain      ! input rotation angle of blade one
   integer, intent(inout) :: iradius      ! number of gridpoints for blade length computed at first call
   real,    intent(inout) :: force(0:ieps,nx,ny,3) ! actuator line force added to force



   real, parameter :: dx=1.0
   real, parameter :: pi=3.1415927410125732
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0
   real, parameter :: eps=2.5              ! smoothing distance in Gaussian kernel 1.6

   real costheta
   real sintheta
   real x0,y0                             ! x-y location of turbine center
   real xb,yb                             ! x-y locatiuon along a blade
   real x,y                               ! x-y locatiuon of a grid point
   real xp,yp                             ! x-y locatiopn of point on blade closest to x-y
   real theta                             ! work blade rotation angle
   real q                                 ! Dynamic pressure
   real ux,vx,wx,dens                     ! local velocity and density variables
   real utheta                            ! tangential velocity from v and w
   real urel2                             ! relative velocity squared
   real, save :: radius                   ! turbine radius in grid cells
   integer i,ia,ib                        ! counter and loop limits in horizontal direction in rotor plane
   integer j,ja,jb                        ! counter and loop limits in vertical direction in rotor plane
   integer k                              ! counter normal to the rotor plane
   integer ic,jc                          ! gridpoint along a blade
   integer iblade                         ! blade counter
   real u0

   real phi                               ! Angle between rotorplane and relative velocity direction
   real cosphi(nrchords)
   real sinphi(nrchords)
   real clift(nrchords)                   ! Lift coefficient for each chord along the blade
   real cdrag(nrchords)                   ! Drag coefficient for each chord along the blade
   real forceL(nrchords)                  ! Lift force along the blade
   real forceD(nrchords)                  ! Drag force along the blade
   real forceT(nrchords)                  ! Tangential (driving force) along the blade
   real forceN(nrchords)                  ! Nornal (drag force) along the blade

   real :: gauss                          ! weigthing function for smeering force
   integer ichord                         ! chord counter
   integer, save :: ifirst=0
   real fL,fD,fN,fT,fTa,fNa

   real, save :: omega                    ! Rotation speed in radians per second
   real, save :: omegand                  ! Nondimensional rotation speed
   real, save :: relma(nrchords),relmb(nrchords)
   character(len=3) tag3
   real newton,R
   real c1,c2,blades,frac,tiploss
   real, save :: g,lambda
   real angattack                         ! computed angle of attack
   real angle                             ! work angle
   integer iave

   ifirst=ifirst+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! I N I T I A L I Z A T I O N   D O N E   A T   F I R S T   C A L L



   if (ifirst == 1) then
! Read lift and drag from foil data files
      call nrelreadfoil()

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
      radius=rotorradius+hubradius     ! radius comes from mod_nrel5Mw
      R=2.0*radius                     ! rotor diameter
      iradius=nint(radius)             ! approximate radius in number of gridcells
      print '(a,f8.2,i4)','Rotor radius=   ',radius,iradius
      print '(a,f8.2)',   'Rotor diamenter=',R
      print '(a,17f8.2)', 'dc  :',dc(1:17)
      print '(a,17f8.2)', 'relm:',relm(1:17)
      print '(a,17f8.2)', 'ra  :',relma(1:17)
      print '(a,17f8.2)', 'chor:',chord(1:17)

! Rotation speed in radians/s
      omega=pi2*turbrpm/60.0
! nondimesonal omega
      omegand=omega*p2l%time

      print *,'Omega  (RPM)                 =',turbrpm
      print *,'Omega  (radians/s)           =',omega
      print *,'Omega  (radians/)            =',omegand
      print *,'Time per revolution          =',pi2/omega
! Tipspeed(Tip speed can be determined from the rotational speed, which is ωR where ω is the rotational
! speed in radians per second and R is the radius of the turbine in meters.)
      print *,'tipspeed R*Omega m/s         =',real(iradius)*p2l%length*omega
! Tipspeed ratio 8 m/s winds (For three blades, a TSR of 6 to 7 is optimal. If it is less, not all
! the available energy is captured; if it is more, the blades move into an area of turbulence from the last
! blade and are not as efficient.

      newton=(p2l%length**2)*p2l%rho*(p2l%vel**2)
      lambda=radius*p2l%length*omega/(uini*p2l%vel)

! Prandtl-Glauert and Shen tip speed model
      if (itiploss == 2) then
         c1=0.125    ! coefficient 1
         c2=21.0     ! coefficient 2
         blades=3.0  ! Number of blades
         g=exp(-c1*(blades*lambda-c2))+0.1
      else
         g=1.0
      endif
   endif

! If given tipspeed ratio is larger than 0.0 we recompute the rotation speed to match the imposed tipspeed ratio
   if (tipspeedratio > 0.0) then
      if (ifirst==1) print '(a)','Recomputing omega based on average  u velocity in rotor plane'
      u0=0.0
      iave=0
      do k=0,360,10
         angle=real(k)*pi2/360.0
         i=ipos+nint(real(iradius)*cos(angle)/2.0)
         j=jpos+nint(real(iradius)*sin(angle)/2.0)
         u0=u0+u(i,j)
         iave=iave+1
      enddo
      u0=u0/real(iave)
      omega=tipspeedratio*u0*p2l%vel/(radius*p2l%length)
      omegand=omega*p2l%time
      lambda=radius*p2l%length*omega/(u0*p2l%vel)
   else
      u0=uini
   endif

   if (ifirst==1) then
      print *
      print '(a)','**************************************'
      print '(a,2f8.2)','TIPSPEED RATIO R*Omega/uini  =',lambda,tipspeedratio
      print '(a,f8.2)','UINI                  (m/s)  =',u0*p2l%vel
      print '(a,f8.2)','omega                 (RPM)  =',omega*60.0/pi2
      print '(a)','**************************************'
      print *
   endif


! Prandtl-Glauert and Shen tip speed model
      if (itiploss == 2) g=exp(-c1*(blades*lambda-c2))+0.1

! The given twist angles (in degrees) refers to zero angle a the tip, and then twists the blade as we
! move towards the hub from chord to chord. We store the twist angle for each chord in twist(:).
! Below we will calculate the angle between the relative velocity and the rotor plane as phi=atan2(ux,utheta).
! In addition there can be a user specified pitchangle.
! The local angle of attach, angattack, is then the difference between phi and the accumulated twist+pitch angles.

! Twist angle at tip
!      angle=atan(1.0/lambda)*360.0/pi2
!      print *,'tip relative angle',angle

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C O M P U T I N G   T H E   B L A D E   F O R C E S
! Set rotation angle
   theta=thetain

! Center point of turbine
   x0=real(ipos)
   y0=real(jpos)

! Loop limits for computing force
   ia=max(1,ipos-iradius-2)
   ib=min(nx,ipos+iradius+2)
   ja=max(1,jpos-iradius-2)
   jb=min(ny,jpos+iradius+2)

! Running through blades
   do iblade=1,3

      costheta=cos(theta)
      sintheta=sin(theta)

! Compute lift and drag force along blade
      do ichord=1,nrchords
         ! Finding gridpoint index closest to chord number ichord to extract velocity
         xb=x0+relm(ichord)*costheta
         yb=y0+relm(ichord)*sintheta
         ic=nint(xb)
         jc=nint(yb)

! Bilinear interpolation, to find u(xb,yb), v(xb,yb), w(xb,yb) in the square  ic=int(xb), ic+1,  jc=int(yb), jc+1
! Non-dimensional velocity components from model
         call bilin(xb,yb,u,v,w,rho,ic,jc,nx,ny,ux,vx,wx,dens)

! Non-dimensional utheta (non-dim omegand and distance relm(ichord)
         utheta =  omegand*relm(ichord) - wx*costheta - vx*sintheta

! Non-dimensional urel**2
         urel2  = ux**2 + utheta**2

! Non-dimensional dynamic pressure
         q= 0.5 * dens * urel2

! Flowangle between relative windspeed and rotor plane
         phi = atan2(ux,utheta)
         sinphi(ichord)=sin(phi)
         cosphi(ichord)=cos(phi)

         angattack=(phi*360.0/pi2-twist(ichord)-pitchangle)
         call nrelliftdrag(clift,cdrag,angattack,ichord)

! Tip loss
         if (itiploss == 0) then
            tiploss=1.0
         else
            blades=3.0
            frac=blades*(radius-relm(ichord))/(2.0*relm(ichord)*sinphi(ichord)+0.00001)
            tiploss=(2.0/pi)*acos(g*exp(-frac))
         endif

! Non-dimensional lift and drag forces per blade element
         forceL(ichord) = q * chord(ichord) * dc(ichord) * clift(ichord) * tiploss
         forceD(ichord) = q * chord(ichord) * dc(ichord) * cdrag(ichord) * tiploss

         forceT(ichord) =forceL(ichord)*sinphi(ichord) - forceD(ichord)*cosphi(ichord)
         forceN(ichord) =forceL(ichord)*cosphi(ichord) + forceD(ichord)*sinphi(ichord)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! D I A G N O S T I C S   F O R   A   B L A D E
         if (ifirst==1) then
            if (iblade==1) then
               if (ichord == 1) then
                  open(10,file='ALMdata.dat')
                  write(10,'(a)')'VARIABLES = "ichord" "dist" "ic" "jc" "ux" "vx" "wx" "utheta" "twist" "phi" "aattack" "urel" &
                                &"chord" "dc" "clift" "cdrag" "tiploss" "forceL(N/m)" "forceD(N/m)" &
                                &"forceT(N/m)" "forceN(N/m)" "forceT()" "forceN()"'
                  write(10,'(a,i4,a)')'ZONE  T="ALMdata" F=Point, I=  ',nrchords,', J=   1, K=1'
               endif

! Lift and drag forces force(LD) per chord in Newton scaled by chord length (DC) to get Newton/meter.
               fL=forceL(ichord)*newton/(dc(ichord)*p2l%length)
               fD=forceD(ichord)*newton/(dc(ichord)*p2l%length)

! Tangential (rotational force) and normal (drag forces) f(TN) in N/m as in zho19a
               fT=fL*sinphi(ichord)  - fD*cosphi(ichord)
               fN=fL*cosphi(ichord)  + fD*sinphi(ichord)

! Non-dimensional tangent and normal  forces using asm20a scaling
               fTa=fT/(p2l%rho*(uini*p2l%vel)**2*radius*p2l%length)
               fNa=fN/(p2l%rho*(uini*p2l%vel)**2*radius*p2l%length)

               write(10,'(i8,f8.4,2i4,50f12.4)')ichord,relm(ichord)/radius,ic,jc,ux,vx,wx,utheta,&
                                                twist(ichord),phi*360.0/pi2,angattack,sqrt(urel2),&
                                                chord(ichord),dc(ichord),&
                                                clift(ichord),cdrag(ichord),&
                                                tiploss, fL, fD, fT, fN, fTa, fNa
               if (ichord == nrchords) close(10)
               !if (ichord == nrchords) stop
            endif
         endif
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing the forces on the grid
!$OMP PARALLEL DO PRIVATE(i,j,ichord,k,y,x,xp,yp,gauss)  &
!$OMP             SHARED(ja,jb,ia,ib,x0,y0,costheta,sintheta,ieps,force,forceN,forceT)
      do j=ja,jb
         y=real(j)
         do i=ia,ib
            x=real(i)
            do ichord=1,nrchords
               xp=x0+relm(ichord)*costheta
               yp=y0+relm(ichord)*sintheta
               do k=0,ieps
                  if (ieps==0) then
                     gauss=(1.0/(pi*eps**2)) * exp(-((x-xp)**2 + (y-yp)**2)/eps**2)
                  else
                     gauss=(1.0/(sqrt(pi**3)*eps**3)) * exp(-((x-xp)**2 + (y-yp)**2 + real(k)**2)/eps**2)
                  endif
                  force(k,i,j,1) = force(k,i,j,1) + forceN(ichord)  * gauss
                  force(k,i,j,2) = force(k,i,j,2) - forceT(ichord)  * sintheta * gauss
                  force(k,i,j,3) = force(k,i,j,3) + forceT(ichord)  * costheta * gauss
               enddo
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO

      theta=theta+rad120
   enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dumping detailed diagnostics for one turbine case
   if ((ifirst <= 0).and.(nturbines==1)) then
      do k=0,ieps
         write(tag3,'(i3.3)')k
         if (ifirst==1) then
            open(21,file='dist'//tag3//'.dat',status='unknown')
            write(21,*)'TITLE = "Dist',tag3,'"'
            write(21,*)'VARIABLES = "i-index" "j-index" "u" "v" "w" "u2"'
            close(21)
         endif

         open(21,file='dist'//tag3//'.dat',position='append')
            write(21,'(3(a,i5),a)')' ZONE  T="k=',k,'" F=BLOCK, I=',nx,', J=',ny,', K=1'
            write(21,'(30I5)')((i,i=1,nx),j=1,ny)
            write(21,'(30I5)')((j,i=1,nx),j=1,ny)
            write(21,'(10(1x,e12.5))')((force(k,i,j,1),i=1,nx),j=1,ny)
            write(21,'(10(1x,e12.5))')((force(k,i,j,2),i=1,nx),j=1,ny)
            write(21,'(10(1x,e12.5))')((force(k,i,j,3),i=1,nx),j=1,ny)
            write(21,'(10(1x,e12.5))')((sqrt(force(k,i,j,2)**2+force(k,i,j,3)**2),i=1,nx),j=1,ny)
         close(21)
      enddo
   endif

end subroutine
end module
