module m_turbines_compute_force
use mod_dimensions
contains
subroutine turbines_compute_force(force,forceN,forceT,thetain,jpos,kpos,iradius)
   use mod_dimensions,  only : ny,nz
   use m_turbines_init, only : ieps
   use mod_nrel5mw,     only : nrchords,relm
   implicit none
   real,    intent(out) :: force(0:ieps,ny,nz,3) ! actuator line force added to force
   real,    intent(in)  :: forceT(nrchords,3)  ! Tangential (driving force) along the blade
   real,    intent(in)  :: forceN(nrchords,3)  ! Nornal (drag force) along the blade
   real,    intent(in)  :: thetain
   integer, intent(in)  :: jpos
   integer, intent(in)  :: kpos
   integer, intent(in)  :: iradius
   real, parameter :: pi=3.1415927410125732
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0
   real, parameter :: eps=2.5

   integer iblade
   integer ichord
   real theta
   real gauss
   integer ja,jb,ka,kb
   real y, z
   real y0,z0
   real yp,zp
   real costheta, sintheta
   integer i,j,k

! Loop limits for computing force
   ja=max(1,jpos-iradius-2)
   jb=min(ny,jpos+iradius+2)
   ka=max(1,kpos-iradius-2)
   kb=min(nz,kpos+iradius+2)
   y0=real(jpos)
   z0=real(kpos)
   theta=thetain


   do iblade=1,3
      costheta=cos(theta)
      sintheta=sin(theta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computing the forces on the grid
!$OMP PARALLEL DO PRIVATE(i,j,k,ichord,y,z,yp,zp,gauss)  &
!$OMP             SHARED(ja,jb,ka,kb,y0,z0,costheta,sintheta,force,forceN,forceT)
      do k=ka,kb
         z=real(k)
         do j=ja,jb
            y=real(j)
            do ichord=1,nrchords
               yp=y0+relm(ichord)*costheta
               zp=z0+relm(ichord)*sintheta
               do i=0,ieps
                  if (ieps==0) then
                     gauss=(1.0/(pi*eps**2)) * exp(-((y-yp)**2 + (z-zp)**2)/eps**2)
                  else
                     gauss=(1.0/(sqrt(pi**3)*eps**3)) * exp(-((y-yp)**2 + (z-zp)**2 + real(i)**2)/eps**2)
                  endif
                  force(i,j,k,1) = force(i,j,k,1) + forceN(ichord,iblade)  * gauss
                  force(i,j,k,2) = force(i,j,k,2) - forceT(ichord,iblade)  * sintheta * gauss
                  force(i,j,k,3) = force(i,j,k,3) + forceT(ichord,iblade)  * costheta * gauss
               enddo
            enddo
         enddo
      enddo
!$OMP END PARALLEL DO

      theta=theta+rad120
   enddo

end subroutine
end module

