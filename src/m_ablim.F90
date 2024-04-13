module m_ablim
contains
subroutine ablim(k,nz,dzfac,ka,kb)
   implicit none
   integer, intent(in)  :: k
   integer, intent(in)  :: nz
   real,    intent(out) :: dzfac
   integer, intent(out) :: ka
   integer, intent(out) :: kb

   if (k==1) then
      ka=1
      kb=k+1
      dzfac=1.0
   elseif (k==nz) then
      ka=k-1
      kb=nz
      dzfac=1.0
   else
      ka=k-1
      kb=k+1
      dzfac=2.0
   endif
end subroutine
end module
