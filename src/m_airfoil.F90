module m_airfoil

contains
subroutine airfoil(blanking)
   use mod_dimensions
   implicit none
   logical, intent(inout)  :: blanking(0:nx+1,0:ny+1,0:nz+1)
   integer, parameter :: num_points = 100
   real :: xcc(num_points)
   real :: ycu(num_points)
   real :: ycl(num_points)
   real :: ycc(num_points)
   real :: ytt(num_points)

! Example NACA 2412 airfoil
! https://en.wikipedia.org/wiki/NACA_airfoil
   real, parameter :: M = 2.0/100.0       ! Curving of chamber (M=0 gives symmetrical chamber)
   real, parameter :: P = 40.0/100.0      ! Location of maximum chamber thickness
   real, parameter :: T = 12.0/100.0      ! Thickness in percent of the chord
   real, parameter :: chord_len =100.0    ! Chord len of the airfoil
   real, parameter :: yscale=150.0        ! Thickness of the airfoil
   real, parameter :: xref=50.0          ! x-starting point of airfoil
   real, parameter :: yref=50.0           ! y-center point of airfoil
   real, parameter :: tilt=-0.10          ! y-center point of airfoil

   integer :: i,j
   real :: theta, yt, yc, dyc_dx, x

   do i = 1, num_points
      x = (i-1)/real(num_points-1)*chord_len

      yt = 5.0*T*(0.2969*sqrt(x/chord_len)&
                - 0.1260*(x/chord_len) &
                - 0.3516*(x/chord_len)**2&
                + 0.2843*(x/chord_len)**3&
                - 0.1015*(x/chord_len)**4)

      if ((x/chord_len) <= P) then
         yc     = (M/P**2)*(2.0*P*(x/chord_len) - (x/chord_len)**2)
         dyc_dx = (2.0*M/P**2) * (P-(x/chord_len))
      else
         yc     = (M/(1.0-P)**2)*( (1.0-2.0*P) + 2.0*P*(x/chord_len) - (x/chord_len)**2)
         dyc_dx = ((2.0*M)/(1.0-P)**2) * (P - (x/chord_len))
      end if

      theta = atan2(dyc_dx, 1.0)

      xcc(i) = xref + x
      ycc(i) = yref + yscale*yc
      ytt(i) = yref + yscale*yt
      ycu(i) = yref + tilt*real(i-1) + yscale*( yc + yt * cos(theta))
      ycl(i) = yref + tilt*real(i-1) + yscale*( yc - yt * cos(theta))
   end do

   ! Print coordinates
   open(10,file='airfoil.dat')
      do i = 1, num_points
         write(10,'(i3,5f10.4)')i,  xcc(i), ytt(i), ycc(i), ycu(i), ycl(i)
      end do
   close(10)

   do j=1,ny
!   do i=nint(xcc(1)),nint(xcc(num_points))
   do i=1,num_points
      if ((nint(ycl(i)) < j) .and. (j <  nint(ycu(i)))) blanking(nint(xref)+i,j,:) = .true.
   enddo
   enddo

!   do j=nint(yref)-25,nint(yref)+25
!      print '(i3,tr2,10001l1)',j,blanking(nint(xcc(1))-10:nint(xcc(num_points))+10,j)
!   enddo

   open(10,file='blanking.dat')
      do j=1,ny
      do i=1,nx
         if (blanking(i,j,1)) write(10,'(2i4,a)')i,j,' 1.0'
      enddo
      enddo
   close(10)
end subroutine
end module


