module m_actuatorline
contains
subroutine actuatorline(forc,nx,ny,ipos,jpos,radius,thetain)
   implicit none
   integer, intent(in)    :: nx           ! dimension of horizontal axis of turbine plane
   integer, intent(in)    :: ny           ! dimension of vertical   axis of turbine plane
   integer, intent(in)    :: ipos         ! horizontal gridpoint for center of turbine
   integer, intent(in)    :: jpos         ! vertical   gridpoint for center of turbine


   real,    intent(in)    :: radius       ! turbine radius
   real,    intent(in)    :: thetain      ! input rotation angle of blade one
   real,    intent(inout) :: forc(nx,ny)  ! actuator line force added to forc


   real, parameter :: dx=1.0
   real, parameter :: pi=3.1415926535
   real, parameter :: pi2=2.0*pi
   real, parameter :: rad120=pi2*120.0/360.0

   real x0,y0                             ! x-y locatiuon of a grid point
   real x1,y1                             ! x-y location of turbine center
   real x2,y2                             ! x-y location of blade tip
   real xp,yp                             ! x-y locatiopn of point on blade closest to x0-y0
   real a,b,t                             ! work variables
   real theta                             ! work rotation angle
   integer iradius                        ! number of gridpoints for blade length
   integer i,ia,ib                        ! counter and loop limits in horizontal direction
   integer j,ja,jb                        ! counter and loop limits in vertical direction
   integer iblade                         ! blade counter

   iradius=nint(radius/dx)

   theta=thetain

! Center point of turbine
   x1=real(ipos)
   y1=real(jpos)

! Loop limits for computing force
   ia=ipos-iradius
   ib=ipos+iradius
   ja=jpos-iradius
   jb=jpos+iradius

!   open(10,file='dist.dat',status='unknown')
!   write(10,*)'TITLE = "Dist"'
!   write(10,*)'VARIABLES = "i-index" "j-index" "Force"'

! Running through blades
   do iblade=1,3

      x2=x1+radius*cos(theta)
      y2=y1+radius*sin(theta)

      a=x2-x1
      b=y2-y1

      do j=ja,jb
         y0=real(j)
         do i=ia,ib
            x0=real(i)

            t=((x0-x1)*a + (y0-y1)*b)/radius**2
            xp=x1+t*a
            yp=y1+t*b

            if ((0.0 <= t).and.(t < 1.0)) then
               forc(i,j) = exp(-((xp-x0)**2 + (yp-y0)**2)/0.5)
            endif

         enddo
      enddo

      theta=theta+rad120
   enddo

!   write(10,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',nx,', J=',ny,', K=1'
!   write(10,'(30I5)')((i,i=1,nx),j=1,ny)
!   write(10,'(30I5)')((j,i=1,nx),j=1,ny)
!   write(10,'(10(1x,e12.5))')((forc(i,j),i=1,nx),j=1,ny)
!   close(10)

end subroutine
end module
