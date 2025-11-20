module m_turbine_rotor_basis
   implicit none
contains
!    Build orthonormal basis for the rotor:
!      e_axis : rotor axis (unit)
!      e1, e2 : in-plane orthonormal vectors
!
!        yaw  : rotation about z-axis
!        tilt : tilt of rotor axis
!--------------------------------------------------------------
#ifdef _CUDA
attributes(host,device) &
#endif
subroutine turbine_rotor_basis(yaw, tilt, e_axis, e1, e2)
   implicit none
   real, intent(in)  :: yaw, tilt
   real, intent(out) :: e_axis(3), e1(3), e2(3)

   e_axis = (/ cos(yaw)*cos(tilt),  sin(yaw)*cos(tilt), -sin(tilt) /)
   e1     = (/ -sin(yaw),           cos(yaw),            0.0       /)
   e2     = (/ -cos(yaw)*sin(tilt), -sin(yaw)*sin(tilt), -cos(tilt) /)

end subroutine
end module


