program LatticeBoltzmann
#ifdef _CUDA
   use cudafor
#endif
   use m_airfoil
   use m_averaging_full
   use m_averaging_sec
   use m_boundarycond
   use m_city
   use m_city2
   use m_collisions
   use m_compute_f
   use m_compute_fneq
   use m_cube
   use m_cylinder
   use m_density
   use m_diag
   use m_drift
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
   use m_velocity
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
   real    :: f(nl,0:nx+1,0:ny+1,0:nz+1)            ! density function
   real    :: feq(nl,0:nx+1,0:ny+1,0:nz+1)          ! Maxwells equilibrium density function
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
#endif

   logical, dimension(:,:,:), allocatable :: lblanking_h ! Host blanking variable
   logical, dimension(:,:,:), allocatable :: lblanking   ! blanking solids grid points
#ifdef _CUDA
   attributes(device) :: lblanking
#endif

   real uvel(nz)
   real, dimension(:), allocatable :: uvel_d        ! vertical u-velocity profile device version
#ifdef _CUDA
   attributes(device) :: uvel_d
#endif

   logical :: lsolids=.false.

! Spatially dependent relaxation time
   real, dimension(:,:,:), allocatable  :: tau      ! relaxation time scale
#ifdef _CUDA
   attributes(device) :: tau
#endif

! Fluid variables
   real    :: u(nx,ny,nz)                           ! x component of fluid velocity
   real    :: v(nx,ny,nz)                           ! y component of fluid velocity
   real    :: w(nx,ny,nz)                           ! z component of fluid velocity
   real    :: rho(nx,ny,nz)                         ! fluid density
#ifdef _CUDA
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: rho
#endif


   real :: elevation(nx,ny)=0.0
   integer i,j,k,l
   integer :: it

   logical, parameter :: debug=.false.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef _CUDA
   print*, "Running GPU version!"
#elif OPEN_MP
   print*, "Running OPEN_MP version!"
#else
   print*, "Running single CPU version!"
#endif

#ifdef D3Q19
   print*, "D3Q19 lattice"
#else
   print*, "D3Q27 lattice"
#endif

#ifdef DOUBLE_PRECISION
   print*, "Double precision code (-r8)"
#else
   print*, "Single precision code"
#endif

   call cpustart()

   allocate(tau(0:nx+1,0:ny+1,0:nz+1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading all input parameters
   call readinfile()

#ifdef _CUDA
   call gpu_meminfo('ini')
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call hermite_polynomials()


   if (nturbines > 0)      call turbines_init()
   if (inflowturbulence)   call inflow_turbulence_init()
#ifdef _CUDA
   call gpu_meminfo('turb')
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      allocate(lblanking(0:nx+1,0:ny+1,0:nz+1))
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k = 0, nz+1
      do j = 0, ny+1
      do i = 0, nx+1
         lblanking(i,j,k) = .false.
      end do
      end do
      end do


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

   if (lsolids) then
      allocate(lblanking_h(0:nx,0:ny,0:nz))
      lblanking_h=lblanking
      do j=1,ny
      do i=1,nx
         do k=1,nz
            if (lblanking_h(i,j,k)) then
                elevation(i,j)=k
            else
                exit
            endif
         enddo
      enddo
      enddo
      call tecfld('elevation',nx,ny,1,elevation)
      deallocate(lblanking_h)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization requires specification of u,v,w, and rho to compute feq
   call uvelshear(uvel)
   allocate(uvel_d(nz))
   uvel_d=uvel

! setting seed to seed.orig if file exist and nt0=0, otherwise generate new seed
   call seedmanagement(nt0)

! setting shapiro factors (no really used)
!   call shfact(nshapiro,sh)

   if (nt0 == 0) then
! Intialization of macro variables
      call inipert(rho,u,v,w,uvel_d)

! Generate trubulence forcing fields
      if (inflowturbulence) call inflow_turbulence_compute(uu,vv,ww,rr,.true.,nrturb)

! Initial diagnostics
      call diag(2,0,rho,u,v,w,lblanking) ! Save initial condidion

! Inititialization with equilibrium distribution from u,v,w, and rho
      call fequil3(feq,rho,u,v,w)
      call boundarycond(feq,uvel_d)
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
      do i=0,nx+1
         do l=1,nl
            f(l,i,j,k) = feq(l,i,j,k)
         enddo
      enddo
      enddo
      enddo

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
      do i=0,nx+1
         !tau(i,j,k) = tauin
         tau(i,j,k) = 3.0*kinevisc + 0.5
      enddo
      enddo
      enddo

   else
! Restart from restart file
      call readrestart(nt0,f,theta,uu,vv,ww,rr)
      call macrovars(rho,u,v,w,f)
! To recover initial tau
      call fequil3(feq,rho,u,v,w)
      call boundarycond(feq,uvel_d)
      call regularization(f, feq, u, v, w)
      if (ivreman == 1) call vreman(f, tau)

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
      do i=0,nx+1
         f(:,i,j,k)=feq(:,i,j,k)+f(:,i,j,k)
      enddo
      enddo
      enddo

   endif
   call diag(1,0,rho,u,v,w,lblanking) ! Save geometry
   call cpufinish(1)

#ifdef _CUDA
   call gpu_meminfo('time stepping')
#endif
   print *,'Start timestepping loop'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation Main Loop
   do it = nt0+1, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) then
         write(*,'(a,i6,a,f10.2,a,3(a,f12.7))',advance='yes')'Iteration:', it,                     &
                                                            ' Time:'  ,real(it)*p2l%time,' s.'
      endif

! start with f,rho,u,v,w


! [feq] = fequil3(rho,u,v,w] (returns equilibrium density)
      call fequil3(feq,rho,u,v,w)

! [f=Rneqf] = regularization[f,feq,u,v,w] (input f is full f and returns reg. non-eq-density)
      call compute_fneq(f,feq)
      if (ihrr == 1) then
         call regularization(f, feq, u, v, w)
      endif

! [u,v,w,turbine_df] = turbines_forcing[rho,u,v,w]
      if (nturbines > 0)      call turbines_forcing(rho,u,v,w)

! [u,v,w,turbulence_df] = turbulenceforcing[rho,u,v,w]
      if (inflowturbulence)   call inflow_turbulence_forcing(rho,u,v,w,turbulence_ampl,it,nrturb)

!! [feq] = fequil3(rho,u,v,w] (returns new equilibrium density using updated forcing velocities)
      if (iforce /= 10 .and. nturbines > 0) then
         call compute_f(f, feq)
         call fequil3(feq,rho,u,v,w)
         call compute_fneq(f, feq)
      endif

! [tau] = vreman[f] [f=Rneqf]
      if (ivreman == 1) call vreman(f, tau)

! [feq=f] = collisions(f,feq,tau)  f=f^eq + (1-1/tau) * R(f^neq)
      call collisions(f,feq,tau)

! [feq=f] = turbines_apply(feq,turbine_df,tau)  f=f+turbine_df
      if (nturbines > 0)      call turbines_apply(feq,turbine_df,tau)

! [feq=f] = inflow_turbulence_apply(feq,turbulence_df,tau)  f=f+turb_df
      if (inflowturbulence)   call inflow_turbulence_apply(feq,turbulence_df)

! Bounce back boundary on fixed walls within the fluid
      if (lsolids) call solids(feq,lblanking)

! General boundary conditions
      call boundarycond(feq,uvel_d)

! Drift of feq returned in f
      call drift(f,feq)


! Compute updated macro variables
      call macrovars(rho,u,v,w,f)

! Diagnostics
      call diag(2,it,rho,u,v,w,lblanking)

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
      if (mod(it,irestart) == 0)            call saverestart(it,f,theta,uu,vv,ww,rr)
      call cpufinish(15)

      if (lmeasurements .and. mod(it,1000)==0) call predicted_measurements(u,v,w,it)

   enddo


   call cpustart()
   call saverestart(it-1,f,theta,uu,vv,ww,rr)
   call cpufinish(16)
   call cpuprint()

   if (allocated(uu))        deallocate(uu)
   if (allocated(vv))        deallocate(vv)
   if (allocated(ww))        deallocate(ww)
   if (allocated(rr))        deallocate(rr)
   if (allocated(lblanking)) deallocate(lblanking)
   if (allocated(uvel_d))    deallocate(uvel_d)

end program LatticeBoltzmann
