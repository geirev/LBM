module m_actuatorline
contains
subroutine actuatorline(forceN,forceT,jpos,kpos,thetain,iradius,u,v,w,rho)
! This routine currently runs on the CPU, as for a single turbine there are only 3 blades time 17 chords, i.e. 51 parallel
! computations. In the future it will make sense to include many turbines (>> 100) in this routine and extend forceT(ichord,iblade)
! to forceT(ichord,iblade,nturbines).
   use mod_dimensions, only : ny,nz
   use m_readinfile, only : p2l,uini,turbrpm,pitchangle,itiploss,tipspeedratio
   use mod_nrel5mw
   use m_nrelreadfoil
   use m_nrelliftdrag
   use m_bilin
   use m_turbines_print_blade

   implicit none
   integer, intent(in)    :: jpos         ! horizontal gridpoint for center of turbine
   integer, intent(in)    :: kpos         ! vertical   gridpoint for center of turbine


   real,    intent(in)    :: u(ny,nz)     ! u velocity through the turbine section
   real,    intent(in)    :: v(ny,nz)     ! v velocity at the turbine section
   real,    intent(in)    :: w(ny,nz)     ! w velocity at the turbine section
   real,    intent(in)    :: rho(ny,nz)   ! density at the turbine section

   real,    intent(in)    :: thetain      ! input rotation angle of blade one
   integer, intent(inout) :: iradius      ! number of gridpoints for blade length computed at first call

   real, intent(out) :: forceT(nrchords,3)  ! Tangential (driving force) along the blade
   real, intent(out) :: forceN(nrchords,3)  ! Nornal (drag force) along the blade


   integer, parameter :: ialmout=2000     ! timestep for dumping ALMdata file
   real, parameter :: pi=3.1415927410125732
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0

   real costheta,sintheta
   real y0,z0                             ! y-z location of turbine center
   real yb,zb                             ! y-z locatiuon along a blade
   integer jc,kc                          ! gridpoint along a blade
   real theta                             ! blade rotation angle
   real dynpres                           ! Dynamic pressure
   real ux,vx,wx,dens                     ! local velocity and density variables
   real utheta                            ! tangential velocity from v and w
   real urel2                             ! relative velocity squared

   integer i,j,k                          ! counters
   integer iblade                         ! blade counter
   integer ichord                         ! chord counter
   real u0                                ! u-velocity at turbine plane

   real phi                               ! Angle between rotorplane and relative velocity direction
   real cosphi(nrchords)
   real sinphi(nrchords)
   real clift(nrchords)                   ! Lift coefficient for each chord along the blade
   real cdrag(nrchords)                   ! Drag coefficient for each chord along the blade
   real forceL(nrchords,3)                ! Lift force along the blade
   real forceD(nrchords,3)                ! Drag force along the blade

! Lift and drag forces force(LD) per chord in Newton scaled by chord length (DC) to get Newton/meter.
!   real fL(nrchords)                      ! for diagnostics
!   real fD(nrchords)                      ! for diagnostics

! Tangential (rotational force) and normal (drag forces) f(TN) in N/m as in zho19a
!   real fN(nrchords)                      ! for diagnostics
!   real fT(nrchords)                      ! for diagnostics

! Non-dimensional tangent and normal  forces using asm20a scaling
!   real fNa(nrchords)                     ! for diagnostics
!   real fTa(nrchords)                     ! for diagnostics

   integer, save :: ifirst=0              ! Counter for calls to routine

   real, save :: radius                   ! turbine radius in grid cells
   real, save :: omega                    ! Rotation speed in radians per second
   real, save :: omegand                  ! Nondimensional rotation speed

   real newton                            ! conversion to Newton

! local tiploss variables
   real c1,c2,blades,frac,tiploss,denom, acos_arg
   real, save :: g                        ! tiploss parameter
   real, save :: lambda
   real angattack                         ! computed angle of attack
   real angle                             ! work angle
   integer iave                           ! counter for computing an average

   ifirst=ifirst+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INITIALIZATION DONE AT FIRST CALL

   if (ifirst == 1) then
! Read lift and drag from foil data files
      call nrelreadfoil()

! Nondimensional turbine parameters
      relm=relm/p2l%length
      dc=dc/p2l%length
      chord=chord/p2l%length
      rotorradius=rotorradius/p2l%length
      hubradius=hubradius/p2l%length

! Blade length in number of gridpoints
      radius=rotorradius+hubradius     ! radius comes from mod_nrel5Mw
      iradius=nint(radius)             ! approximate radius in number of gridcells
      print '(a,f8.2,i4)','Rotor radius=   ',radius,iradius
      print '(a,f8.2)',   'Rotor diamenter=',2.0*radius
!      print '(a,50f8.2)', 'dc  :',dc(1:nrchords)
!      print '(a,50f8.2)', 'relm:',relm(1:nrchords)
!      print '(a,50f8.2)', 'chor:',chord(1:nrchords)

! Rotation speed in radians/s
      omega=pi2*turbrpm/60.0
! nondimesonal omega
      omegand=omega*p2l%time

!      print *,'Omega  (RPM)                 =',turbrpm
!      print *,'Omega  (radians/s)           =',omega
!      print *,'Omega  (radians/)            =',omegand
      !print *,'Time per revolution          =',pi2/omega
! Tipspeed(Tip speed can be determined from the rotational speed, which is ωR where ω is the rotational
! speed in radians per second and R is the radius of the turbine in meters.)
      print *,'tipspeed R*Omega m/s         =',real(iradius)*p2l%length*omega
! Tipspeed ratio 8 m/s winds (For three blades, a TSR of 6 to 7 is optimal. If it is less, not all
! the available energy is captured; if it is more, the blades move into an area of turbulence from the last
! blade and are not as efficient.

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

! If given tipspeed ratio  from infile.in is larger than 0.0 we recompute the rotation speed to match the imposed tipspeed ratio
! based on the average velocity sampled on the circle at half a rotor radius.
   if (tipspeedratio > 0.0) then
      if (ifirst==1) print '(a)','Recomputing omega and tipspeed ratio based on average  u velocity in rotor plane'
      u0=0.0
      iave=0
      do i=0,360,10
         angle=real(k)*pi2/360.0
         j=jpos+nint(real(iradius)*cos(angle)/2.0)
         k=kpos+nint(real(iradius)*sin(angle)/2.0)
         u0=u0+u(j,k)
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
      print '(a)','********************************************************'
      if (tipspeedratio > 0.0) then
         print '(a)','using prescribed tip-speed-ratio and recomputed omega'
      else
         print '(a)','tip-speed-ratio computed from given omega and u-inflow'
      endif
      print '(a,f8.2,i4)','Rotor radius          (m,i)  =',radius*p2l%length,iradius
      print '(a,f8.2)',   'Rotor diamenter       (m)    =',2.0*radius*p2l%length
      print '(a,f8.2)',   'TIPSPEED RATIO R*Omega/uini  =',lambda
      print '(a,f8.2)',   'UINI                  (m/s)  =',u0*p2l%vel
      print '(a,f8.2)',   'omega                 (RPM)  =',omega*60.0/pi2
      print '(a,f8.2)',   'Time per revolution   (s)    =',pi2/omega
      print '(a)','********************************************************'
      print *
   endif


! Prandtl-Glauert and Shen tip speed model (updated with new lambda)
      if (itiploss == 2) g=exp(-c1*(blades*lambda-c2))+0.1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! C O M P U T I N G   T H E   B L A D E   F O R C E S
! The given twist angles (in degrees) refers to zero angle a the tip, and then twists the blade as we
! move towards the hub from chord to chord. We store the twist angle for each chord in twist(:).
! Below we will calculate the angle between the relative velocity and the rotor plane as phi=atan2(ux,utheta).
! In addition there can be a user specified pitchangle.
! The local angle of attach, angattack, is then the difference between phi and the accumulated twist+pitch angles.

! Twist angle at tip
!      angle=atan(1.0/lambda)*360.0/pi2
!      print *,'tip relative angle',angle

! Set rotation angle

! Center point of turbine
   y0=real(jpos)
   z0=real(kpos)

! Running through blades
   theta=thetain
   do iblade=1,3

      costheta=cos(theta)
      sintheta=sin(theta)

! Compute lift and drag force along blade
      do ichord=1,nrchords
         ! Finding pivot gridpoint index closest to chord number ichord to extract velocity
         yb=y0+relm(ichord)*costheta
         zb=z0+relm(ichord)*sintheta
         jc=int(yb)
         kc=int(zb)

! Bilinear interpolation, to find u(yb,zb), v(yb,zb), w(yb,zb) in the square  jc=int(yb), jc+1,  kc=int(zb), kc+1
! Non-dimensional velocity components from model
         call bilin(yb,zb,u,v,w,rho,jc,kc,ny,nz,ux,vx,wx,dens)

! Non-dimensional utheta (non-dim omegand and distance relm(ichord)
         utheta =  omegand*relm(ichord) - wx*costheta - vx*sintheta

! Non-dimensional urel**2
         urel2  = ux**2 + utheta**2

! Non-dimensional dynamic pressure
         dynpres= 0.5 * dens * urel2

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
            denom = 2.0*relm(ichord)*max(abs(sinphi(ichord)), 1.0e-6)
            frac = blades*(radius - relm(ichord)) / denom
            acos_arg = g * exp(-frac)
            if (acos_arg >  1.0) acos_arg =  1.0
            if (acos_arg < -1.0) acos_arg = -1.0
            tiploss = (2.0/pi) * acos(acos_arg)
         endif

! Non-dimensional lift and drag forces per blade element
         forceL(ichord,iblade) = dynpres * chord(ichord) * dc(ichord) * clift(ichord) * tiploss
         forceD(ichord,iblade) = dynpres * chord(ichord) * dc(ichord) * cdrag(ichord) * tiploss

         forceT(ichord,iblade) =forceL(ichord,iblade)*sinphi(ichord) - forceD(ichord,iblade)*cosphi(ichord)
         forceN(ichord,iblade) =forceL(ichord,iblade)*cosphi(ichord) + forceD(ichord,iblade)*sinphi(ichord)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Diagnostocs
!         if (iblade == 1 .and. ifirst == ialmout) then
!         ! Lift and drag forces force(LD) per chord in Newton scaled by chord length (DC) to get Newton/meter.
!            newton=(p2l%length**2)*p2l%rho*(p2l%vel**2)
!            fL(ichord)=forceL(ichord,iblade)*newton/(dc(ichord)*p2l%length)
!            fD(ichord)=forceD(ichord,iblade)*newton/(dc(ichord)*p2l%length)
!
!         ! Tangential (rotational force) and normal (drag forces) f(TN) in N/m as in zho19a
!            fT(ichord)=fL(ichord)*sinphi(ichord)  - fD(ichord)*cosphi(ichord)
!            fN(ichord)=fL(ichord)*cosphi(ichord)  + fD(ichord)*sinphi(ichord)
!
!         ! Non-dimensional tangent and normal  forces using asm20a scaling
!            fTa(ichord)=fT(ichord)/(p2l%rho*(uini*p2l%vel)**2*radius*p2l%length)
!            fNa(ichord)=fN(ichord)/(p2l%rho*(uini*p2l%vel)**2*radius*p2l%length)
!         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      enddo
      theta=theta+rad120
   enddo

!   if (ifirst == ialmout) then
! To recover results from Asmuth 2020 (fig 4):
!    Run with constant RPM of 8.95 (tipspeedratio set to 0.0 meaning constant RPM),
!    use inflow velocity of 8.00 m/s,
!    and save result after a significant spinup.
!      call turbines_print_blade(clift,cdrag,fL,fD,fN,fT,fTa,fNa,radius)
!   endif

end subroutine
end module


