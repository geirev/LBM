module m_nrelliftdrag

contains
subroutine nrelliftdrag(cl,cd,angle,ichord)
   use mod_nrel5mw
   implicit none
   real, intent(inout)  :: cl(nrchords) ! Lift coefficients
   real, intent(inout)  :: cd(nrchords) ! Drag coefficients
   real, intent(in)     :: angle        ! angle between vrel and rotor plane
   integer, intent(in)  :: ichord

   real slope
   integer i,j,k

   j=ichord     ! computing lift and drag for chord nr ichord
   k=nfoil(j)   ! k refers to the data file for chord(ichord)

   do i=1,nrang(k)
      if ( foil(i,k)%deg < angle .and. angle <= foil(i+1,k)%deg ) then
         slope=(foil(i+1,k)%cl-foil(i,k)%cl)/(foil(i+1,k)%deg-foil(i,k)%deg)
         cl(j)=foil(i,k)%cl + slope*(angle - foil(i,k)%deg)

         slope=(foil(i+1,k)%cd-foil(i,k)%cd)/(foil(i+1,k)%deg-foil(i,k)%deg)
         cd(j)=foil(i,k)%cd + slope*(angle - foil(i,k)%deg)

         print '(2i4,3(a,3f10.4))',j,i , ' angle: ', foil(i,k)%deg , angle, foil(i+1,k)%deg, &
                       ' -- ',  foil(i,k)%cl, cl(j), foil(i+1,k)%cl,&
                       ' -- ',  foil(i,k)%cd, cd(j), foil(i+1,k)%cd
         exit
      endif
   enddo


end subroutine

end module
