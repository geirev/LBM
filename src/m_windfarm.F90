module m_windfarm
contains
subroutine windfarm(blanking)
   use mod_dimensions
   logical, intent(inout) :: blanking(nx,ny,nz)
   integer i,j,k

! No-slip at ground
   blanking(:,:,1)=.true.

end subroutine
end module

