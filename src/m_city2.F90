module m_city2
contains
subroutine city2(blanking)
   use mod_dimensions, only : nx, nyg, nz
   implicit none
   logical, intent(inout) :: blanking(0:nx+1,0:nyg+1,0:nz+1)
   integer, parameter :: nrb=20
   integer :: ipos(nrb)
   integer :: jpos(nrb)
   integer i,j,kpos,irad,ib

   ipos(1)=60; jpos(1)=18
   ipos(2)=64; jpos(2)=jpos(1)+16
   ipos(3)=61; jpos(3)=jpos(2)+10
   ipos(4)=59; jpos(4)=jpos(3)+18
   ipos(5)=62; jpos(5)=jpos(4)+16

   ipos(6)=78; jpos(6)=19
   ipos(7)=78; jpos(7)=jpos(6)+12
   ipos(8)=79; jpos(8)=jpos(7)+14
   ipos(9)=76; jpos(9)=jpos(8)+14
   ipos(10)=74; jpos(10)=jpos(9)+15

   ipos(11)=98; jpos(11)=21
   ipos(12)=93; jpos(12)=jpos(11)+12
   ipos(13)=91; jpos(13)=jpos(12)+14
   ipos(14)=99; jpos(14)=jpos(13)+17
   ipos(15)=90; jpos(15)=jpos(14)+18

   ipos(16)=111; jpos(16)=17
   ipos(17)=114; jpos(17)=jpos(16)+12
   ipos(18)=109; jpos(18)=jpos(17)+14
   ipos(19)=110; jpos(19)=jpos(18)+14
   ipos(20)=116; jpos(20)=jpos(19)+18

   jpos(:)=jpos(:)+12

! Building row one
   irad=4
   kpos=100
   do ib=1,nrb
      do j=0,nyg+1
      do i=0,nx+1
         if ( (abs(i-ipos(ib)) < irad) .and. abs(j-jpos(ib)) < irad ) then
            blanking(i,j,1:min(kpos,nz)) = .true.
         endif
      enddo
      enddo
   enddo

end subroutine
end module
