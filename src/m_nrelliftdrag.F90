module m_nrelliftdrag

contains
subroutine nrelliftdrag(cl,cd,angle,ichord)
   use mod_nrel5mw
   implicit none
   real, intent(inout)  :: cl(nrchords) ! Lift coefficients
   real, intent(inout)  :: cd(nrchords) ! Drag coefficients
   real, intent(in)     :: angle        ! angle between vrel and rotor plane
   integer, intent(in)  :: ichord       ! ichord number

   real slope
   integer i,k

   k=nfoil(ichord)   ! k refers to the data file for chord(ichord)

   do i=1,nrang(k)-1
      if ( foil(i,k)%deg < angle .and. angle <= foil(i+1,k)%deg ) then
         slope=(foil(i+1,k)%cl-foil(i,k)%cl)/(foil(i+1,k)%deg-foil(i,k)%deg)
         cl(ichord)=foil(i,k)%cl + slope*(angle - foil(i,k)%deg)

         slope=(foil(i+1,k)%cd-foil(i,k)%cd)/(foil(i+1,k)%deg-foil(i,k)%deg)
         cd(ichord)=foil(i,k)%cd + slope*(angle - foil(i,k)%deg)

!         print '(3i4,3(a,3f10.4))',k,ichord,i , ' angle: ', foil(i,k)%deg , angle, foil(i+1,k)%deg, &
!                       ' -- ',  foil(i,k)%cl, cl(ichord), foil(i+1,k)%cl,&
!                       ' -- ',  foil(i,k)%cd, cd(ichord), foil(i+1,k)%cd
         exit
      endif
   enddo


end subroutine

end module
