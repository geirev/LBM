module m_closedbnd
contains
subroutine closedbnd(blanking)
   use mod_dimensions, only : nx,ny,nz
   use m_readinfile,   only : ibndbb,jbndbb,kbndbb
   logical, intent(inout) :: blanking(nx,ny,nz)

! No-slip bounce-back boundaries
   if (ibndbb == 11 .or. ibndbb == 33) blanking( 1, :, :)=.true.
   if (ibndbb == 22 .or. ibndbb == 33) blanking(nx, :, :)=.true.

   if (jbndbb == 11 .or. jbndbb == 33) blanking( :, 1, :)=.true.
   if (jbndbb == 22 .or. jbndbb == 33) blanking( :,ny, :)=.true.

   if (kbndbb == 11 .or. kbndbb == 33) blanking( :, :, 1)=.true.
   if (kbndbb == 22 .or. kbndbb == 33) blanking( :, :,nz)=.true.

end subroutine
end module

