program LatticeBoltzmann
#ifdef _CUDA
   use cudafor
#endif
   use m_airfoil
   use m_testing
   use m_postcoll
   use m_printdefines
   use m_averaging_full
   use m_averaging_sec
   use m_fillghosts
   use m_boundarycond
   use m_city
   use m_city2
   use m_collisions
   use m_compute_f
   use m_compute_fneq
   use m_cube
   use m_cylinder
   use m_diag
   use m_drift
   use m_dump_elevation
   use m_fequil3
   use m_gpu_meminfo
   use m_inflow_turbulence_apply
   use m_inflow_turbulence_compute
   use m_inflow_turbulence_forcing
   use m_inflow_turbulence_init
   use m_inipert
   use m_macrovars
   use m_predicted_measurements
   use m_readinfile
   use m_readrestart
   use m_regularization
   use m_rhotest
   use m_saverestart
   use m_seedmanagement
   use m_set_random_seed3
   use m_solids
   use m_sphere
   use m_tecfld
   use m_turbines_apply
   use m_turbines_forcing
   use m_turbines_init
   use m_uvelshear
   use m_vreman
   use m_wtime
   use mod_D3Q27setup
   use mod_dimensions
   use mod_shapiro
   use, intrinsic :: omp_lib
   implicit none

!   integer, parameter :: nshapiro=4
!   real sh(0:nshapiro)

! Main variables
   real, target ::     fA(nl,0:nx+1,0:ny+1,0:nz+1) ! density function
   real, target ::     fB(nl,0:nx+1,0:ny+1,0:nz+1) ! equilibrium density function
   real         ::       tau(0:nx+1,0:ny+1,0:nz+1) ! relaxation time scale
   logical      :: lblanking(0:nx+1,0:ny+1,0:nz+1) ! blanking solids grid points
   real         ::         u(nx,ny,nz)             ! x component of fluid velocity
   real         ::         v(nx,ny,nz)             ! y component of fluid velocity
   real         ::         w(nx,ny,nz)             ! z component of fluid velocity
   real         ::       rho(nx,ny,nz)             ! fluid density
   real         ::      uvel(nz)                   ! vertical u-velocity profile on device

   real, pointer:: f1(:,:,:,:), f2(:,:,:,:), tmp(:,:,:,:)

#ifdef _CUDA
   attributes(device) :: f1
   attributes(device) :: f2
   attributes(device) :: tmp
   attributes(device) :: fA
   attributes(device) :: fB
   attributes(device) :: tau
   attributes(device) :: lblanking
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: rho
   attributes(device) :: uvel
#endif

   real,    dimension(:),       allocatable :: uvel_h      ! vertical u-velocity profile on host
   integer :: it
   logical :: lsolids=.false.


   call printdefines()
   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading all input parameters
   call readinfile()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some initialization
   call hermite_polynomials()
   if (nturbines > 0)      call turbines_init()
   if (inflowturbulence)   call inflow_turbulence_init()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Print GPU memory use
#ifdef _CUDA
   call gpu_meminfo('ini')
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Including solid bodies like cylinder and urban city as blanking of grid cells
   lblanking = .false.

   select case(trim(experiment))
   case('city')
      call city(lsolids,lblanking)
   case('city2')
      call city2(lsolids,lblanking)
   case('cylinder')
      call cylinder(lsolids,lblanking)
   case('airfoil')
!      call airfoil(lsolids,lblanking)
       stop 'needs fix airfoil routine for gpu'
   end select
   if (lsolids) call dump_elevation(lblanking)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save grid geometry and blanking
   call diag(1,0,rho,u,v,w,lblanking)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setting seed to seed.orig if file exist and nt0=0, otherwise generate new seed
   call seedmanagement(nt0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setting shapiro factors (no really used)
!   call shfact(nshapiro,sh)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization requires specification of u,v,w, and rho to compute feq

! Inflow velocity shear stored in uvel as read from file in uvelshare
   allocate(uvel_h(nz))
   call uvelshear(uvel_h)
   uvel=uvel_h
   deallocate(uvel_h)

   if (nt0 == 0) then
! Intialization of macro variables
      call inipert(rho,u,v,w,uvel)

! Initial diagnostics
      call diag(itecout,0,rho,u,v,w,lblanking)

! Inititialization with equilibrium distribution from u,v,w, and rho
      call fequil3(fB,rho,u,v,w)
      call fillghosts(fB)
      call boundarycond(fB,fA,uvel)
      fA=fB

! Generate turbulence forcing fields
      if (inflowturbulence) call inflow_turbulence_compute(uu,vv,ww,rr,.true.,nrturb)
   else
! Restart from restart file
      call readrestart(nt0,fA,theta,uu,vv,ww,rr)
      call macrovars(rho,u,v,w,fA)
! To ensure we have values in the boundary points first time we call boundarycond.
      fB=fA
   endif
   tau = 3.0*kinevisc + 0.5 ! static and constant tau also for ghost zones in case they are used

   call cpufinish(1)

#ifdef _CUDA
   call gpu_meminfo('time stepping')
#endif


! initialize pointers
   f1 => fA
   f2 => fB



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   print *,'Start timestepping loop'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do it = nt0+1, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) then
         write(*,'(a,i6,a,f10.2,a,3(a,f12.7))')'Iteration:', it,' Time:'  ,real(it)*p2l%time,' s.'
      endif

! start time step with f1,rho,u,v,w given

! [turbine_df = turbines_forcing(rho,u,v,w)]
      if (nturbines > 0)      call turbines_forcing(rho,u,v,w,it)

! [turbulence_df = turbulenceforcing(rho,u,v,w)]
      if (inflowturbulence)   call inflow_turbulence_forcing(rho,u,v,w,turbulence_ampl,it,nrturb)

! [f1,tau = post collision(f1,rho,u,v,w] (returns post collision density and tau for forcing)
      call postcoll(f1, tau, rho,u,v,w)

! [f1 = f1 + turbine_df]
      if (nturbines > 0)      call turbines_apply(f1,turbine_df,tau)

! [f1 = f1 + turbulence_df]
      if (inflowturbulence)   call inflow_turbulence_apply(f1,turbulence_df)

! Bounce back boundary on fixed walls within the fluid
      if (lsolids) call solids(f1,lblanking)

! [f1 and f2 updated with boundary conditions]
      call boundarycond(f1,f2,uvel)


! Drift of f1 returned in f2
      call drift(f2,f1)

! Swap such that f2 becomes f1 in next time step
      tmp   => f1
      f1  => f2
      f2 => tmp

! Compute updated macro variables
      call macrovars(rho,u,v,w,f1)

! Diagnostics
      call diag(itecout,it,rho,u,v,w,lblanking)

      call cpustart()
! Averaging for diagnostics similar to Asmuth paper for wind turbines
      if (laveturb .and. nturbines > 0) then
         if (avestart < it .and. it < avesave) call averaging_sec(u,v,w,.false.,iradius)
         if (it == avesave)                    call averaging_sec(u,v,w,.true.,iradius)
      endif

! Averaging for diagnostics
      if (laveraging) then
         if (avestart < it .and. it < avesave) call averaging_full(u,v,w,rho,lblanking,.false.)
         if (it == avesave)                    call averaging_full(u,v,w,rho,lblanking,.true.)
      endif

! Updating input turbulence matrix
      if (mod(it, nrturb) == 0 .and. it > 1 .and. inflowturbulence .and. ibnd==1) then
         call inflow_turbulence_compute(uu,vv,ww,rr,.false.,nrturb)
      endif

! Save restart file
      if (mod(it,irestart) == 0)            call saverestart(it,f1,theta,uu,vv,ww,rr)
      call cpufinish(15)

      if (lmeasurements .and. mod(it,1000)==0) call predicted_measurements(u,v,w,it)

   enddo


   call cpustart()
   call saverestart(it-1,f1,theta,uu,vv,ww,rr)
   call cpufinish(16)
   call cpuprint()

   if (allocated(uu       ))  deallocate(uu            )
   if (allocated(vv       ))  deallocate(vv            )
   if (allocated(ww       ))  deallocate(ww            )
   if (allocated(rr       ))  deallocate(rr            )

   call testing(it,f1,f2)

end program LatticeBoltzmann
