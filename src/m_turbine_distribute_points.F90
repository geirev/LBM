module m_turbine_distribute_points
   implicit none
contains

! Build the global list of actuator sample points (points_global) from turbine
! hub locations and blade geometry.
!
! Runtime switch per turbine, via turbines_in(it)%imodel:
!   imodel = 0 : actuator line (ALM)   - nblades discrete points, rotating with theta
!   imodel = 1 : rotating actuator disk (ADM-R) - nazim static points per radius,
!                each carrying force_scale = nblades/nazim so the annulus receives
!                the same total force as nblades individual blade passes would.
subroutine turbine_distribute_points(turbines_in, points_global)

   use mod_turbines,    only : turbine_t, point_t, pi2
   use mod_turbine_def, only : nrchords, relm, dc, chord, twist
   use m_turbine_extend_array
   use m_turbine_rotor_basis

   implicit none

   type(turbine_t),              intent(in)  :: turbines_in(:)
   type(point_t),   allocatable, intent(out) :: points_global(:)

   integer :: it, ib, ic
   integer :: np
   integer :: nsamples
   real    :: e_axis(3), e1(3), e2(3), e_rot(3)
   real    :: theta, yaw, tilt
   real    :: fscale
   type(point_t) :: pt

   np = 0
   allocate(points_global(np))

   do it = 1, size(turbines_in)

      yaw  = turbines_in(it)%yaw
      tilt = turbines_in(it)%tilt

      call turbine_rotor_basis(yaw, tilt, e_axis, e1, e2)

      if (turbines_in(it)%imodel == 0) then

         if (turbines_in(it)%nblades <= 0) then
            write(*,*) 'ERROR: nblades must be positive for turbine ', it
            stop
         end if

         nsamples = turbines_in(it)%nblades
         fscale   = 1.0

      else if (turbines_in(it)%imodel == 1) then

         if (turbines_in(it)%nazim <= 0) then
            write(*,*) 'ERROR: nazim must be positive for ADM-R turbine ', it
            stop
         end if

         nsamples = turbines_in(it)%nazim
         fscale   = real(turbines_in(it)%nblades) / &
                    real(turbines_in(it)%nazim)

      else
         write(*,*) 'ERROR: invalid turbine model ', turbines_in(it)%imodel
         stop
      end if

      do ib = 1, nsamples

         if (turbines_in(it)%imodel == 0) then
            ! Blade azimuth advances with rotor rotation.
            theta = turbines_in(it)%theta + &
                    real(ib-1)*pi2/real(nsamples)
         else
            ! Fixed ring of sample points, independent of instantaneous theta.
            theta = real(ib-1)*pi2/real(nsamples)
         end if

         do ic = 1, nrchords

            pt%iturb  = it
            pt%iblade = ib
            pt%ichord = ic

            e_rot(:) = cos(theta)*e1(:) + sin(theta)*e2(:)

            pt%xg = turbines_in(it)%xhub + relm(ic)*e_rot(1)
            pt%yg = turbines_in(it)%yhub + relm(ic)*e_rot(2)
            pt%zg = turbines_in(it)%zhub + relm(ic)*e_rot(3)

            pt%yaw         = yaw
            pt%tilt        = tilt
            pt%theta       = theta

            pt%relm        = relm(ic)
            pt%dc          = dc(ic)
            pt%chord       = chord(ic)
            pt%twist       = twist(ic)
            pt%foil        = ic

            pt%pitch       = turbines_in(it)%pitchangle
            pt%omegand     = turbines_in(it)%omegand
            pt%force_scale = fscale

            call turbine_extend_array(points_global, np+1)
            np = np + 1
            points_global(np) = pt

         end do
      end do
   end do

end subroutine turbine_distribute_points

end module m_turbine_distribute_points
