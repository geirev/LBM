!=======================================================================
! Turbine time-series diagnostics
!
! Writes one ASCII time-series file for each turbine:
!
!      output/turbine_001.dat
!      output/turbine_002.dat
!      ...
!
! One line is appended every time turbine_diagnostics_timeseries()
! is called.
!
! Columns:
!
!   1   timestep       model timestep                           [-]
!   2   time           physical simulation time                [s]
!   3   ulocal         filtered local upstream wind speed      [m/s]
!   4   rpm            rotor speed                             [rev/min]
!   5   TSR            actual tip-speed ratio                  [-]
!   6   Ct             thrust coefficient                      [-]
!   7   Cq             torque coefficient                      [-]
!   8   Cp             power coefficient                       [-]
!   9   T              rotor thrust                            [kN]
!  10   Q              rotor torque                            [kN m]
!  11   P              mechanical rotor power                 [MW]
!
! No header is written to the output files; each line contains values only.
!
! Physical conversion factors:
!
!      force_conv  = rho_c L_c^4/T_c^2                  [N/LB-force]
!      torque_conv = rho_c L_c^5/T_c^2                  [N m/LB-torque]
!      power_conv  = rho_c L_c^5/T_c^3                  [W/LB-power]
!
! The global actuator-point forces have already been combined using
! MPI_Allreduce before this routine is called. Consequently only MPI
! rank 0 writes the files.
!=======================================================================
module m_turbine_diagnostics_timeseries

   implicit none

   logical, save :: files_initialized = .false.

contains

subroutine turbine_diagnostics_timeseries(turbines_in,points_global, &
                                          Fvec_global,np,itimestep)

#ifdef MPI
   use mpi
#endif

   use mod_turbines, only : turbine_t,point_t, &
                            uavg_f,vavg_f,wavg_f, &
                            windfilter_initialized

   use mod_turbine_def, only : radiusnd
   use m_readinfile, only : nturbines,p2l,uini,powerloss
   use m_turbine_rotor_basis

#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif

   implicit none

   type(turbine_t), intent(in) :: turbines_in(:)
   type(point_t),   intent(in) :: points_global(:)
   real,            intent(in) :: Fvec_global(3,np)
   integer,         intent(in) :: np
   integer,         intent(in) :: itimestep

#ifndef MPI
   integer, parameter :: mpi_rank = 0
#endif

   integer :: n,p
   integer :: iu,ios

   real, parameter :: pi      = acos(-1.0)
   real, parameter :: pi2     = 2.0*pi
   real, parameter :: rho_ref = 1.0

   real :: e_axis(3),e1(3),e2(3)
   real :: e_tan(3)
   real :: rvec(3),fpoint(3),torque_vec(3)

   real :: thrust             ! [LB force]
   real :: torque             ! [LB torque]
   real :: tq                 ! [LB torque]
   real :: power              ! [LB power]

   real :: radius             ! [LU]
   real :: area               ! [LU^2]
   real :: ulocal             ! [LU/TS]

   real :: omega              ! [rad/s]
   real :: rpm                ! [rev/min]
   real :: lambda_actual      ! [-]

   real :: ct,cq,cp           ! [-]

   real :: force_conv         ! [N/LB-force]
   real :: torque_conv        ! [N m/LB-torque]
   real :: power_conv         ! [W/LB-power]

   real :: thrust_phys        ! [N]
   real :: torque_phys        ! [N m]
   real :: power_phys         ! [W]
   real :: ulocal_phys        ! [m/s]
   real :: time_phys          ! [s]


   character(len=256) :: filename
   logical ex


!-----------------------------------------------------------------------
! Only rank 0 writes because Fvec_global already contains the result
! of the MPI_Allreduce over all MPI tiles.
!-----------------------------------------------------------------------
   if (mpi_rank /= 0) return


!-----------------------------------------------------------------------
! Physical conversion factors.
!-----------------------------------------------------------------------
   force_conv  = p2l%rho*p2l%length**4/p2l%time**2
   torque_conv = p2l%rho*p2l%length**5/p2l%time**2
   power_conv  = p2l%rho*p2l%length**5/p2l%time**3


!-----------------------------------------------------------------------
! Physical simulation time.
!
! This assumes timestep 1 corresponds to t=0, consistent with
!
!      t = real(itimestep-1)*p2l%time
!
! used elsewhere in the model.
!-----------------------------------------------------------------------
   time_phys = real(itimestep-1)*p2l%time


!-----------------------------------------------------------------------
! Create output directory and initialize output files on first call.
!
! status='replace' ensures a new simulation does not silently append
! data to files left from a previous run.
!-----------------------------------------------------------------------
   if (.not.files_initialized) then

      call execute_command_line('mkdir -p output')

      do n = 1,nturbines

         write(filename,'("output/turbine_",a,".dat")') trim(turbines_in(n)%name)
         inquire(file=trim(filename),exist=ex)

         if ((.not.ex) .or. (ex.and. itimestep == 1)) then
            open(newunit=iu,file=trim(filename),  status='replace',action='write',iostat=ios)
            if (ios /= 0) then
               write(*,*) 'ERROR opening turbine diagnostic file: ', trim(filename)
               stop
            endif
            write(iu,'(2A)')'# ',trim(turbines_in(n)%name)
            write(iu,'(A)',advance='no')'# timestep'
            write(iu,'(A)',advance='no')'            time[s]'
            write(iu,'(A)',advance='no')'             U[m/s]'
            write(iu,'(A)',advance='no')'       RPM[rev/min]'
            write(iu,'(A)',advance='no')'             TSR[-]'
            write(iu,'(A)',advance='no')'         Pitch[deg]'
            write(iu,'(A)',advance='no')'              Ct[-]'
            write(iu,'(A)',advance='no')'              Cq[-]'
            write(iu,'(A)',advance='no')'              Cp[-]'
            write(iu,'(A)',advance='no')'              T[kN]'
            write(iu,'(A)',advance='no')'             Q[kNm]'
            write(iu,'(A)',advance='no')'          Paero[MW]'
            write(iu,'(A)',advance='no')'          Pelec[MW]'
            write(iu,*)
            close(iu)
         endif
      enddo

      files_initialized = .true.

   endif


!=======================================================================
! Turbine loop
!=======================================================================
   do n = 1,nturbines


!-----------------------------------------------------------------------
! Local upstream wind must have been initialized before meaningful
! aerodynamic coefficients can be calculated.
!-----------------------------------------------------------------------
      if (.not.windfilter_initialized(n)) cycle


!-----------------------------------------------------------------------
! Rotor basis.
!-----------------------------------------------------------------------
      call turbine_rotor_basis(turbines_in(n)%yaw, &
                               turbines_in(n)%tilt, &
                               e_axis,e1,e2)


!-----------------------------------------------------------------------
! Integrate thrust and torque from all actuator points belonging to
! this turbine.
!-----------------------------------------------------------------------
      thrust = 0.0
      torque = 0.0

      do p = 1,np

         if (points_global(p)%iturb /= n) cycle

         fpoint(:) = Fvec_global(:,p)

         ! Radius vector [LU].
         rvec(1) = points_global(p)%xg - turbines_in(n)%xhub
         rvec(2) = points_global(p)%yg - turbines_in(n)%yhub
         rvec(3) = points_global(p)%zg - turbines_in(n)%zhub

         ! Tangential unit vector [-].
         e_tan(:) = -sin(points_global(p)%theta)*e1(:) + &
                      cos(points_global(p)%theta)*e2(:)

         ! Axial thrust [LB force].
         thrust = thrust + dot_product(fpoint,e_axis)

         ! Rotor torque Q = (r x F) . e_axis [LB torque].
         torque_vec(1) = rvec(2)*fpoint(3) - rvec(3)*fpoint(2)
         torque_vec(2) = rvec(3)*fpoint(1) - rvec(1)*fpoint(3)
         torque_vec(3) = rvec(1)*fpoint(2) - rvec(2)*fpoint(1)

         tq = dot_product(torque_vec,e_axis)

         torque = torque + tq

      enddo


!-----------------------------------------------------------------------
! Local filtered upstream wind speed.
!-----------------------------------------------------------------------
      ulocal = sqrt(uavg_f(n)**2 + &
                    vavg_f(n)**2 + &
                    wavg_f(n)**2)

      ulocal_phys = ulocal*p2l%vel


!-----------------------------------------------------------------------
! Rotor geometry.
!-----------------------------------------------------------------------
      area   = pi*radiusnd*radiusnd


!-----------------------------------------------------------------------
! Rotor angular velocity.
!-----------------------------------------------------------------------
      omega = turbines_in(n)%omegand/p2l%time
      rpm   = omega*60.0/pi2


!-----------------------------------------------------------------------
! Aerodynamic performance quantities.
!-----------------------------------------------------------------------
      if (ulocal > 1.0e-12) then

         ! Tip-speed ratio [-].
         lambda_actual = turbines_in(n)%omegand*radiusnd/ulocal

         ! Mechanical rotor power [LB power].
         power = torque*turbines_in(n)%omegand

         ! Thrust coefficient [-].
         ct = thrust / &
              (0.5*rho_ref*ulocal**2*area)

         ! Torque coefficient [-].
         cq = torque / &
              (0.5*rho_ref*ulocal**2*area*radiusnd)

         ! Power coefficient [-].
         cp = power / &
              (0.5*rho_ref*ulocal**3*area)

      else

         lambda_actual = 0.0
         power         = 0.0
         ct            = 0.0
         cq            = 0.0
         cp            = 0.0

      endif


!-----------------------------------------------------------------------
! Convert loads to physical SI units.
!-----------------------------------------------------------------------
      thrust_phys = thrust*force_conv       ! [N]
      torque_phys = torque*torque_conv      ! [N m]
      power_phys  = power*power_conv        ! [W]


!-----------------------------------------------------------------------
! Append one ASCII line containing VALUES ONLY:
!
! timestep time[s] U[m/s] RPM TSR Ct Cq Cp T[kN] Q[kN m] P[MW]
!-----------------------------------------------------------------------
      write(filename,'("output/turbine_",a,".dat")') trim(turbines_in(n)%name)

      open(newunit=iu,file=trim(filename), &
           status='old',position='append',action='write',iostat=ios)

      if (ios /= 0) then
         write(*,*) 'ERROR opening turbine diagnostic file: ', &
                    trim(filename)
         stop
      endif

      write(iu,'(I10,1X,12(ES18.10,1X))') &
           itimestep, &
           time_phys, &
           ulocal_phys, &
           rpm, &
           lambda_actual, &
           turbines_in(n)%pitchangle, &
           ct, &
           cq, &
           cp, &
           thrust_phys/1.0e3, &
           torque_phys/1.0e3, &
           power_phys/1.0e6,  &
           power_phys*powerloss/1.0e6

      close(iu)

   enddo

end subroutine turbine_diagnostics_timeseries

end module m_turbine_diagnostics_timeseries
