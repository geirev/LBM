module m_city
contains
subroutine city(blanking)
   use mod_dimensions, only : nx, nyg, nz
   implicit none
   logical, intent(inout) :: blanking(0:nx+1,0:nyg+1,0:nz+1)
   integer i,j,ipos,jpos,kpos,irad,jshift


   ipos=60
   kpos=50
   irad=4

   do ipos=60,110,4*irad
   kpos=kpos-10
   if (jshift > 0) then
      jshift=-jshift-1
   else
      jshift=-jshift+1
   endif
   do jpos=9,90,4*irad
      do j=0,nyg+1
      do i=0,nx+1
         if ( (abs(i-ipos) < irad) .and. abs(j-jpos) < irad ) then
            blanking(i,j,0:min(kpos,nz+1)) = .true.
         endif
      enddo
      enddo
   enddo
   enddo

end subroutine
end module

