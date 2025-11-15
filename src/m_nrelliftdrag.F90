module m_nrelliftdrag
contains

#ifdef _CUDA
attributes(host,device) &
#endif
subroutine nrelliftdrag(cl, cd, angle, k)
   use mod_nrel5mw
   implicit none

   real,    intent(out) :: cl, cd
   real,    intent(in)  :: angle
   integer, intent(in)  :: k

   real :: slope
   integer :: i

   cl = 0.0
   cd = 0.0

   do i = 1, nrang(k)-1
      if (foil(i,k)%deg < angle .and. angle <= foil(i+1,k)%deg) then

         slope = (foil(i+1,k)%cl - foil(i,k)%cl) / &
                 (foil(i+1,k)%deg - foil(i,k)%deg)
         cl = foil(i,k)%cl + slope*(angle - foil(i,k)%deg)

         slope = (foil(i+1,k)%cd - foil(i,k)%cd) / &
                 (foil(i+1,k)%deg - foil(i,k)%deg)
         cd = foil(i,k)%cd + slope*(angle - foil(i,k)%deg)

         return
      end if
   end do

end subroutine nrelliftdrag
end module m_nrelliftdrag
