module m_city
contains
subroutine city(lsolids,blanking)
   use mod_dimensions
   implicit none
   logical, intent(out)   :: lsolids
   logical, intent(inout) :: blanking(0:nx+1,0:ny+1,0:nz+1)
   integer i,j,ipos,jpos,kpos,irad,jshift
#ifdef _CUDA
   attributes(device) :: blanking
#endif
   lsolids=.true.

! Building row one
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
      do j=0,ny+1
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
