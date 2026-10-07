module m_liftdrag
contains


!=======================================================================
! Lift/drag interpolation for GPU/device calculations.
!
! Uses the device copies nrang and foil.
!=======================================================================
#ifdef _CUDA
attributes(device) &
#endif
subroutine liftdrag(cl, cd, angle, k)

   use mod_turbine_def, only : nrang,foil

   implicit none

   real,    intent(out) :: cl, cd
   real,    intent(in)  :: angle
   integer, intent(in)  :: k

   real :: slope
   integer :: i

   cl = 0.0
   cd = 0.0

   do i = 1, nrang(k)-1

      if (foil(i,k)%deg < angle .and. &
          angle <= foil(i+1,k)%deg) then

         slope = (foil(i+1,k)%cl - foil(i,k)%cl) / &
                 (foil(i+1,k)%deg - foil(i,k)%deg)

         cl = foil(i,k)%cl + slope*(angle - foil(i,k)%deg)

         slope = (foil(i+1,k)%cd - foil(i,k)%cd) / &
                 (foil(i+1,k)%deg - foil(i,k)%deg)

         cd = foil(i,k)%cd + slope*(angle - foil(i,k)%deg)

         return

      endif

   enddo

end subroutine liftdrag


!=======================================================================
! Lift/drag interpolation for CPU/host calculations.
!
! Uses the host copies nrang_h and foil_h.
! This routine is used by the near-wake smearing correction, which is
! presently evaluated on the CPU.
!=======================================================================
subroutine liftdrag_host(cl, cd, angle, k)

   use mod_turbine_def, only : nrang_h,foil_h

   implicit none

   real,    intent(out) :: cl, cd
   real,    intent(in)  :: angle
   integer, intent(in)  :: k

   real :: slope
   integer :: i

   cl = 0.0
   cd = 0.0

   do i = 1, nrang_h(k)-1

      if (foil_h(i,k)%deg < angle .and. &
          angle <= foil_h(i+1,k)%deg) then

         slope = (foil_h(i+1,k)%cl - foil_h(i,k)%cl) / &
                 (foil_h(i+1,k)%deg - foil_h(i,k)%deg)

         cl = foil_h(i,k)%cl + slope*(angle - foil_h(i,k)%deg)

         slope = (foil_h(i+1,k)%cd - foil_h(i,k)%cd) / &
                 (foil_h(i+1,k)%deg - foil_h(i,k)%deg)

         cd = foil_h(i,k)%cd + slope*(angle - foil_h(i,k)%deg)

         return

      endif

   enddo

end subroutine liftdrag_host


end module m_liftdrag
