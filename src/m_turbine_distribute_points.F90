module m_turbine_distribute_points
   implicit none
contains

! Build the global list of actuator sample points (points_global) from turbine
! hub locations and blade geometry.
!
! Runtime switch per turbine, via:
!   actuator_model = 0 : actuator line (ALM) - nblades discrete azimuthal positions,
!                        rotating with theta.
!   actuator_model = 1 : rotating actuator disk (ADM-R) - nazim fixed azimuthal
!                        positions per radius, each carrying
!                        force_scale = nblades/nazim so that the annulus receives
!                        the same total force as nblades individual blade passes.
subroutine turbine_distribute_points(turbines_in, points_global)

   use mod_turbines,    only : turbine_t, point_t, pi2
   use mod_turbine_def, only : nrchords, relm, dc, chord, twist, radiusnd, nblades
   use m_turbine_rotor_basis
   use m_readinfile, only : actuator_model,nazim

   implicit none

   type(turbine_t),              intent(in)  :: turbines_in(:)
   type(point_t),   allocatable, intent(out) :: points_global(:)

   integer :: it, ib, ic
   integer :: np, p
   integer :: nsamples

   real :: e_axis(3), e1(3), e2(3), e_rot(3)
   real :: theta, yaw, tilt
   real :: fscale

   type(point_t) :: pt


!-----------------------------------------------------------------------
! Count total number of actuator points.
!-----------------------------------------------------------------------
   np = 0

   do it = 1,size(turbines_in)

      select case(actuator_model)

      case(0)
         if (nblades <= 0) then
            write(*,*) 'ERROR: nblades must be positive for turbine ',it
            error stop
         endif

         nsamples = nblades

      case(1)
         if (nazim <= 0) then
            write(*,*) 'ERROR: nazim must be positive for ADM-R turbine ',it
            error stop
         endif

         nsamples = nazim

      case default
         write(*,*) 'ERROR: invalid turbine model for turbine ',it, &
                    ': ',actuator_model
         error stop

      end select

      np = np + nsamples*nrchords

   enddo


!-----------------------------------------------------------------------
! Allocate the complete actuator-point array once.
!-----------------------------------------------------------------------
   allocate(points_global(np))

   p = 0


!-----------------------------------------------------------------------
! Construct actuator points.
!-----------------------------------------------------------------------
   do it = 1,size(turbines_in)

! Define rotor plane
      yaw  = turbines_in(it)%yaw
      tilt = turbines_in(it)%tilt

      call turbine_rotor_basis(yaw,tilt,e_axis,e1,e2)


! ALM or ADM-R model
      select case(actuator_model)

      case(0)
         nsamples = nblades
         fscale   = 1.0

      case(1)
         nsamples = nazim
         fscale   = real(nblades)/real(nazim)

      end select


! Loop over blade or azimuthal positions
      do ib = 1,nsamples

         if (actuator_model == 0) then

            ! Blade azimuth advances with rotor rotation.
            theta = turbines_in(it)%theta + &
                    real(ib-1)*pi2/real(nsamples)

         else

            ! Fixed ring of ADM-R sampling points.
            theta = real(ib-1)*pi2/real(nsamples)

         endif

         e_rot(:) = cos(theta)*e1(:) + sin(theta)*e2(:)


! Loop over radial blade elements
         do ic = 1,nrchords

            p = p + 1

            pt%iturb  = it
            pt%iblade = ib
            pt%ichord = ic

            pt%xg = turbines_in(it)%xhub + relm(ic)*e_rot(1)
            pt%yg = turbines_in(it)%yhub + relm(ic)*e_rot(2)
            pt%zg = turbines_in(it)%zhub + relm(ic)*e_rot(3)

            pt%yaw   = yaw
            pt%tilt  = tilt
            pt%theta = theta

            pt%relm  = relm(ic)
            pt%dc    = dc(ic)
            pt%chord = chord(ic)
            pt%twist = twist(ic)
            pt%foil  = ic

            pt%pitch       = turbines_in(it)%pitchangle
            pt%omegand     = turbines_in(it)%omegand
            pt%force_scale = fscale

            points_global(p) = pt

         enddo
      enddo
   enddo


!-----------------------------------------------------------------------
! Consistency check.
!-----------------------------------------------------------------------
   if (p /= np) then
      write(*,*) 'ERROR turbine_distribute_points: p /= np'
      write(*,*) 'p  = ',p
      write(*,*) 'np = ',np
      error stop
   endif

end subroutine turbine_distribute_points

end module m_turbine_distribute_points
