module m_turbine_compute_bladeforce
contains

!--------------------------------------------------------------
!  subroutine turbine_compute_blade_force
!
!  PURPOSE:
!    Given axial and tangential relative velocities (u_ax,
!    u_tan_rel) in the rotor plane, plus cl/cd, compute the
!    3D force vector on the blade section in global coordinates.
!
!  ARGUMENTS:
!    point      : point_t containing geometry + orientation
!    u_ax       : axial component of relative velocity
!    u_tan_rel  : tangential component (Ω r - tangential flow)
!    dens       : local fluid density
!    cl, cd     : lift & drag coefficients
!    Fvec(3)    : resulting force vector (global x,y,z)
!--------------------------------------------------------------
#ifdef _CUDA
attributes(host,device) &
#endif
subroutine turbine_compute_bladeforce(Fvec, point, u_ax, u_tan_rel, dens, cl, cd)
   use mod_turbines, only : point_t
   use m_turbine_rotor_basis
   implicit none
   type(point_t), intent(in) :: point
   real, intent(in)          :: u_ax, u_tan_rel
   real, intent(in)          :: dens
   real, intent(in)          :: cl, cd
   real, intent(out)         :: Fvec(3)

   real :: e_axis(3), e1(3), e2(3)
   real :: e_rot(3), e_tan(3)
   real :: speed2, area, L, D
   real :: phi, sinphi, cosphi
   real :: tiploss

   ! Build rotor coordinate system
   call turbine_rotor_basis(point%yaw, point%tilt, e_axis, e1, e2)

   ! Radial and tangential directions
   e_rot =  cos(point%theta)*e1 + sin(point%theta)*e2
   e_tan = -sin(point%theta)*e1 + cos(point%theta)*e2

   ! Relative speed squared in rotor plane
   speed2 = u_ax*u_ax + u_tan_rel*u_tan_rel

   ! Chord area
   area   = point%chord * point%dc

   ! Tip-loss factor (placeholder)
   tiploss = 1.0

   L = 0.5 * dens * speed2 * cl * area * tiploss
   D = 0.5 * dens * speed2 * cd * area * tiploss

   ! Flow angle
   phi    = atan2(u_ax, u_tan_rel)
   sinphi = sin(phi)
   cosphi = cos(phi)

   ! Lift and drag decomposition in (e_axis, e_tan) plane
   Fvec = L * (cosphi * e_axis + sinphi * e_tan) - &
          D * (sinphi * e_axis - cosphi * e_tan)
end subroutine
end module
