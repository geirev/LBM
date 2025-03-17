program testbnd

   implicit none

   real, parameter       :: cs2=1/3.0
   real, parameter       :: cs4=1/9.0
   real, parameter       :: cs6=1/27.0

   integer cxs(-1:1,-1:1,-1:1)
   integer cys(-1:1,-1:1,-1:1)
   integer czs(-1:1,-1:1,-1:1)
   integer i,j,k
   real weigths(-1:1,-1:1,-1:1)

   do k=-1,1
   do j=-1,1
   do i=-1,1
      cxs(i,j,k)=i
      cys(i,j,k)=j
      czs(i,j,k)=k
      if (abs(i)+abs(j)+abs(k) == 0) then
         weigths(i,j,k) = 8.0/27.0
      elseif (abs(i)+abs(j)+abs(k) == 1) then
         weigths(i,j,k) = 2.0/27.0
      elseif (abs(i)+abs(j)+abs(k) == 2) then
         weigths(i,j,k) = 1.0/54.0
      elseif (abs(i)+abs(j)+abs(k) == 3) then
         weigths(i,j,k) = 1.0/216.0
      else
         print *,'problem X'
      endif

   enddo
   enddo
   enddo

   print '(a,f10.4)','8.0/27.0  :',8.0/27.0
   print '(a,f10.4)','2.0/27.0  :',2.0/27.0
   print '(a,f10.4)','1.0/54.0  :',1.0/54.0
   print '(a,f10.4)','1.0/216.0 :',1.0/216.0

   write(*,'(a,27i3)')'cxs:',cxs
   write(*,'(a,27i3)')'cys:',cys
   write(*,'(a,27i3)')'czs:',czs
   write(*,'(a)')'weigths:'
   write(*,'(9f10.4)')weigths

end program

! Define velocity vectors
!                         0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6  7  8  9  0  1  2  3  4  5  6
!   integer :: cxs(1:nl) = [0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1,-1, 1, 0, 0,-1, 1, 0, 0,-1, 1,-1, 1, 1,-1,-1, 1]
!   integer :: cys(1:nl) = [0, 0, 0, 1,-1, 0, 0, 1,-1,-1, 1, 0, 0, 1,-1, 0, 0,-1, 1, 1,-1,-1, 1, 1,-1, 1,-1]
!   integer :: czs(1:nl) = [0, 0, 0, 0, 0,-1, 1, 0, 0, 0, 0,-1, 1, 1,-1, 1,-1, 1,-1, 1,-1,-1, 1,-1, 1,-1, 1]
!
!  ! Define weights
!   real :: weights(1:nl) = [8./27.,                                                                     & ! 0-center
!        2./27., 2./27., 2./27., 2./27., 2./27., 2./27.,                                                 & ! 1-6
!        1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., 1./54., & ! 7-18
!        1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216., 1./216.]                           ! 19-26


