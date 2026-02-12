module m_city3
contains
subroutine city3(blanking)
   use mod_dimensions, only : nx, nyg, nz
   implicit none
   logical, intent(inout) :: blanking(0:nx+1,0:nyg+1,0:nz+1)
   integer ioff
   integer joff

   ioff=10
   joff=0

   blanking(ioff+ 3: ioff+ 4, joff+  3:  joff+ 4, 1:3)=.true.
   blanking(ioff+ 3: ioff+ 4, joff+ 11:  joff+12, 1:5)=.true.
   blanking(ioff+ 3: ioff+ 4, joff+ 19:  joff+20, 1:3)=.true.
   blanking(ioff+ 3: ioff+ 4, joff+ 27:  joff+28, 1:1)=.true.
   blanking(ioff+19: ioff+20, joff+  3:  joff+ 4, 1:5)=.true.
   blanking(ioff+19: ioff+20, joff+ 11:  joff+12, 1:2)=.true.
   blanking(ioff+19: ioff+20, joff+ 19:  joff+20, 1:6)=.true.
   blanking(ioff+19: ioff+20, joff+ 27:  joff+28, 1:3)=.true.
   blanking(ioff+11: ioff+12, joff+  7:  joff+ 8, 1:5)=.true.
   blanking(ioff+11: ioff+12, joff+ 15:  joff+16, 1:3)=.true.
   blanking(ioff+11: ioff+12, joff+ 23:  joff+24, 1:2)=.true.
   blanking(ioff+11: ioff+12, joff+  1:  joff+ 1, 1:3)=.true.
   blanking(ioff+11: ioff+12, joff+ 31:  joff+31, 1:3)=.true.
   blanking(ioff+27: ioff+28, joff+  7:  joff+ 8, 1:3)=.true.
   blanking(ioff+27: ioff+28, joff+ 15:  joff+16, 1:2)=.true.
   blanking(ioff+27: ioff+28, joff+ 23:  joff+24, 1:5)=.true.
   blanking(ioff+27: ioff+28, joff+  1:  joff+ 1, 1:3)=.true.
   blanking(ioff+27: ioff+28, joff+ 31:  joff+31, 1:3)=.true.

end subroutine
end module
