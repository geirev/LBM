module m_windfarm
contains
subroutine windfarm(blanking,kbnd)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer, intent(in)    :: kbnd

! No-slip at ground
   if (kbnd /= 0 ) then
      blanking(:,:,1)=.true.
   endif

end subroutine
end module

