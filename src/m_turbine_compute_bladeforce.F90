module m_turbine_compute_bladeforce
contains

#ifdef _CUDA
attributes(host,device) &
#endif
subroutine turbine_compute_bladeforce(Fvec, point, u_ax, u_tan_rel, dens, cl, cd, rtip)

   use mod_turbines, only : point_t
   use mod_turbine_def, only : nblades
   use m_turbine_rotor_basis

   implicit none

   type(point_t), intent(in) :: point
   real, intent(in)          :: u_ax, u_tan_rel
   real, intent(in)          :: dens
   real, intent(in)          :: rtip
   real, intent(in)          :: cl, cd
   real, intent(out)         :: Fvec(3)

   real :: e_axis(3), e1(3), e2(3)
   real :: e_tan(3)
   real :: speed2, area, L, D
   real :: phi, sinphi, cosphi

   real :: ftip, arg
   real :: r
   real :: exp_arg

   real, parameter :: pi = acos(-1.0)
   real, parameter :: epsloss = 1.0e-6

   ! Build rotor coordinate system
   call turbine_rotor_basis(point%yaw, point%tilt, e_axis, e1, e2)

   ! Tangential direction
   e_tan = -sin(point%theta)*e1 + cos(point%theta)*e2

   ! Relative speed squared
   speed2 = u_ax*u_ax + u_tan_rel*u_tan_rel

   ! Blade-element area
   area = point%chord * point%dc

   ! Flow angle
   phi    = atan2(u_ax, u_tan_rel)
   sinphi = sin(phi)
   cosphi = cos(phi)

   !------------------------------------------------------------
   ! Prandtl tip-loss correction
   !
   ! Ftip = 2/pi acos[exp(-f)]
   !
   ! f = B/2 * (R-r)/(r |sin(phi)|)
   !------------------------------------------------------------
   r    = point%relm

   arg = 0.5*real(nblades) * max(rtip-r,0.0) / &
         max(r*abs(sinphi),epsloss)

   exp_arg = min(1.0,max(0.0,exp(-arg)))

   ftip = (2.0/pi)*acos(exp_arg)

   ftip = max(0.0,min(1.0,ftip))

   ! Aerodynamic forces
   L = 0.5 * dens * speed2 * cl * area * ftip
   D = 0.5 * dens * speed2 * cd * area * ftip

   ! Lift and drag decomposition in (e_axis,e_tan) plane
   Fvec = L * (cosphi*e_axis + sinphi*e_tan) + &
          D * (sinphi*e_axis - cosphi*e_tan)

end subroutine turbine_compute_bladeforce

end module m_turbine_compute_bladeforce
