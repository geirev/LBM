module m_turbines_print_blade
! Dumping detailed diagnostics for one turbine blade
contains
subroutine turbines_print_blade(clift,cdrag,fL,fD,fN,fT,fTa,fNa,radius)
   use mod_dimensions,  only : ny,nz
   use m_turbines_init, only : ieps
   use mod_nrel5mw,     only : nrchords, relm
   implicit none
   real, intent(in) :: radius

   real, intent(in) :: clift(nrchords)
   real, intent(in) :: cdrag(nrchords)

! Lift and drag forces force(LD) per chord in Newton scaled by chord length (DC) to get Newton/meter.
   real, intent(in) :: fL(nrchords)
   real, intent(in) :: fD(nrchords)

! Tangential (rotational force) and normal (drag forces) f(TN) in N/m as in zho19a
   real, intent(in) :: fN(nrchords)
   real, intent(in) :: fT(nrchords)

! Non-dimensional tangent and normal  forces using asm20a scaling
   real, intent(in) :: fNa(nrchords)
   real, intent(in) :: fTa(nrchords)

   integer :: ichord

   open(10,file='ALMdata.dat')
      write(10,'(a)')'VARIABLES = "ichord" "dist" "clift" "cdrag" "forceL(N/m)" "forceD(N/m)" &
                                 &"forceT(N/m)" "forceN(N/m)" "forceT()" "forceN()"'
      write(10,'(a,i4,a)')'ZONE  T="ALMdata" F=Point, I=  ',nrchords,', J=   1, K=1'

      do ichord=1,nrchords
         write(10,'(i8,f8.4,8f12.4)')ichord,relm(ichord)/radius,clift(ichord),cdrag(ichord),&
                                          fL(ichord), fD(ichord), fT(ichord), fN(ichord), fTa(ichord), fNa(ichord)
      enddo

   close(10)
end subroutine
end module
