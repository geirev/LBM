module m_turbine_distribute_points
   implicit none
contains
! Build the global list of actuator sample points (points_global) from turbine hub locations and blade geometry.
subroutine turbine_distribute_points(turbines_in, points_global)
   use mod_turbines, only : turbine_t, point_t, pi2
   use m_turbine_extend_array
   use m_turbine_rotor_basis
#ifdef MPI
   use m_mpi_decomp_init, only : j_start, j_end, mpi_rank
#endif
#ifdef _CUDA
   use cudafor
#endif
   implicit none
   type(turbine_t),              intent(in)  :: turbines_in(:)
   type(point_t),   allocatable, intent(out) :: points_global(:)

   integer :: it, ib, ic
   integer :: np
   real    :: e_axis(3), e1(3), e2(3), e_rot(3)
   real    :: theta, yaw, tilt
   type(point_t) :: pt

   np = 0
   allocate(points_global(np))

   do it = 1, size(turbines_in)
      yaw  = turbines_in(it)%yaw
      tilt = turbines_in(it)%tilt

      call turbine_rotor_basis(yaw, tilt, e_axis, e1, e2)

      do ib = 1, turbines_in(it)%nblades
         theta = turbines_in(it)%theta + real(ib-1)*pi2/real(turbines_in(it)%nblades)

         do ic = 1, turbines_in(it)%nchords
            pt%iturb  = it
            pt%iblade = ib
            pt%ichord = ic

            e_rot(:) = cos(theta)*e1(:) + sin(theta)*e2(:)

            pt%xg = turbines_in(it)%xhub + turbines_in(it)%relm(ic) * e_rot(1)
            pt%yg = turbines_in(it)%yhub + turbines_in(it)%relm(ic) * e_rot(2)
            pt%zg = turbines_in(it)%zhub + turbines_in(it)%relm(ic) * e_rot(3)

            pt%yaw     = yaw
            pt%tilt    = tilt
            pt%theta   = theta
            pt%relm    = turbines_in(it)%relm(ic)
            pt%dc      = turbines_in(it)%dc(ic)
            pt%chord   = turbines_in(it)%chord(ic)
            pt%twist   = turbines_in(it)%twist(ic)
            pt%pitch   = turbines_in(it)%pitchangle
            pt%foil    = turbines_in(it)%nfoil(ic)
            pt%omegand = turbines_in(it)%omegand

            call turbine_extend_array(points_global, np+1)
            np = np + 1
            points_global(np) = pt
         end do
      end do
   end do
end subroutine turbine_distribute_points
end module

