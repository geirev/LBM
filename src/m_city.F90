module m_city
contains
subroutine city(lblanking)
   use mod_dimensions
   use m_tecfld
   logical, intent(inout) :: lblanking(0:nx+1,0:ny+1,0:nz+1)
   integer i,j,ipos,jpos,kpos,irad
   real elevation(nx,ny)

! Building row one
   elevation=0.0
   ipos=60
   kpos=50
   irad=4
   do ipos=60,110,4*irad
   kpos=kpos-10
   do jpos=9,90,4*irad
      do j=0,ny+1
      do i=0,nx+1
         if ( (abs(i-ipos) < irad) .and. abs(j-jpos) < irad ) then
            lblanking(i,j,0:kpos) = .true.
            elevation(i,j)=real(kpos)
         endif
      enddo
      enddo
   enddo
   enddo


   call tecfld('elevation',nx,ny,1,elevation)


end subroutine
end module

