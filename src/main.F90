program LatticeBoltzmann
#ifdef _CUDA
   use cudafor
#endif
   use m_mechanical_ablvisc
   use m_abl_initialize
   use m_advection
   use m_heatflux
   use m_buoyancy_forcing
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
   use m_inflow_turbulence_update
   use m_inflow_turbulence_forcing
   use m_inflow_turbulence_init
   use mod_turbines
   use m_turbine_initialize
   use m_turbine_forcing
   use m_forcings_apply
   use m_inipert
   use m_macrovars
   use m_predicted_measurements
   use m_readinfile
   use m_readrestart
   use m_regularization
   use m_rhotest
   use m_saverestart
   use m_save_uvw
   use m_seedmanagement
   use m_solids
   use m_solid_objects_init
   use m_sphere
   use m_tecfld
   use m_uvelshear
   use m_vreman
   use m_wtime
   use mod_D3Q27setup
   use mod_dimensions
   use mod_shapiro
   use, intrinsic :: omp_lib

#ifdef MPI
   use mpi
   use m_mpi_decomp_init
   use m_mpi_halo_exchange_j
   use m_mpi_decomp_finalize
#endif

   implicit none

!   integer, parameter :: nshapiro=4
!   real sh(0:nshapiro)

! Main variables
   real, target ::     fA(nl,0:nx+1,0:ny+1,0:nz+1) ! density function
   real, target ::     fB(nl,0:nx+1,0:ny+1,0:nz+1) ! equilibrium density function
   !real         ::     fH(nl,0:nx+1,0:ny+1,0:nz+1) ! equilibrium density function
   real         ::       tau(0:nx+1,0:ny+1,0:nz+1) ! relaxation time scale
   logical      :: lblanking(0:nx+1,0:ny+1,0:nz+1) ! blanking solids grid points
   real         ::         u(0:nx+1,0:ny+1,0:nz+1)             ! x component of fluid velocity
   real         ::         v(0:nx+1,0:ny+1,0:nz+1)             ! y component of fluid velocity
   real         ::         w(0:nx+1,0:ny+1,0:nz+1)             ! z component of fluid velocity
   real         ::       rho(0:nx+1,0:ny+1,0:nz+1)             ! fluid density
   real         ::      uvel(nz)                   ! vertical u-velocity profile on device
   real         ::         u_h(0:nx+1,0:ny+1,0:nz+1)             ! x component of fluid velocity
   real         ::         v_h(0:nx+1,0:ny+1,0:nz+1)             ! y component of fluid velocity
   real         ::         w_h(0:nx+1,0:ny+1,0:nz+1)             ! z component of fluid velocity
   real         ::       rho_h(0:nx+1,0:ny+1,0:nz+1)             ! fluid density

   real, target, allocatable  :: tracerA(:,:,:,:)
   real, target, allocatable  :: tracerB(:,:,:,:)

   real, target, allocatable  :: pottempA(:,:,:)
   real, target, allocatable  :: pottempB(:,:,:)

   real,         allocatable  :: external_forcing(:,:,:,:)

   real, pointer:: f1(:,:,:,:), f2(:,:,:,:), tmp(:,:,:,:)
   real, pointer:: t1(:,:,:,:), t2(:,:,:,:), tr_tmp(:,:,:,:)
   real, pointer:: p1(:,:,:), p2(:,:,:), pt_tmp(:,:,:)

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

   attributes(device) :: tracerA,tracerB,t1,t2,tr_tmp
   attributes(device) :: pottempA,pottempB,p1,p2,pt_tmp
   attributes(device) :: external_forcing
#endif

   real,    dimension(:),       allocatable :: uvel_h      ! vertical u-velocity profile on host
   integer :: it,ir,i
   logical :: lsolids=.false.

#ifdef MPI
   logical periodic_j_bc
#endif


   call printdefines()
   call cpustart()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading all input parameters
   call readinfile()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some additional allocations

   if (ntracer > 0) then
      allocate(tracerA(ntracer,0:nx+1,0:ny+1,0:nz+1))
      allocate(tracerB(ntracer,0:nx+1,0:ny+1,0:nz+1))
   endif

   if (iablvisc == 2) then
      allocate(pottempA(0:nx+1,0:ny+1,0:nz+1))
      allocate(pottempB(0:nx+1,0:ny+1,0:nz+1))
   endif

   allocate(external_forcing(3,0:nx+1,0:ny+1,0:nz+1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialize MPI and domain decomposition
   ir=0
#ifdef MPI
   periodic_j_bc = (jbnd == 0)
   call mpi_decomp_init(periodic_j_bc)
   if (mpi_rank == 0) print *, 'MPI ranks=', mpi_nprocs, ' nyg=', nyg, ' ny=', ny
   ir=mpi_rank
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Some initialization
   call hermite_polynomials()
   if (nturbines > 0)      call turbine_initialize()
   if (inflowturbulence)   call inflow_turbulence_init()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Including solid bodies like cylinder and urban city as blanking of grid cells
   call solid_objects_init(lblanking, lsolids, trim(experiment), ir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Save grid geometry and blanking
   call diag(1,0,rho,u,v,w,pottempA,tracerA,lblanking)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setting seed to seed.orig if file exist and nt0=0, otherwise generate new seed
   call seedmanagement(nt0,ir)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! setting shapiro factors (not used)
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
      call inipert(rho,u,v,w,uvel,ir)

! tracer initialization
      if (ntracer >= 1) then
            tracerA(1,0:nx+1,0:ny+1,0:nz+1)=1.0
            tracerA(1,100:110,110:120,0:nz+1)=2.0
            tracerA(1,150:160,50:60,0:nz+1)=2.0
      endif

! Initialization if mechanical or full atmospheric boundary layer
      call mechanical_ablvisc(ir)
      if (iablvisc == 2) then
         call abl_initialize(pottempA,ir)
      endif

#ifdef MPI
      if (ntracer > 0) then
         call mpi_halo_exchange_j(tracerA,ntracer)
      endif
      if (iablvisc == 2) then
         call mpi_halo_exchange_j(pottempA,1)
      endif
#endif

! Initial diagnostics
      call diag(itecout,0,rho,u,v,w,pottempA,tracerA,lblanking)

      print *,'finsihed diag',ir
! Inititialization with equilibrium distribution from u,v,w, and rho
      call fequil3(fB,rho,u,v,w)
      call fillghosts(fB)
      fA=fB
      call boundarycond(fB,fA,uvel,tracerA,pottempA)
      call boundarycond(fA,fB,uvel,tracerA,pottempA)
      if (ntracer > 0) tracerB=tracerA
      if (iablvisc == 2) pottempB=pottempA

! Generate turbulence forcing fields
      if (inflowturbulence) call inflow_turbulence_update(uu,vv,ww,rr,nrturb,.false.)
   else
! Restart from restart file
      call readrestart(nt0,fA,turbines(1:nturbines)%theta,uu,vv,ww,rr,pottempA,tracerA)
      call macrovars(rho,u,v,w,fA)
! To ensure we have values in the boundary points first time we call boundarycond.
      fB=fA
      if (ntracer > 0) tracerB=tracerA
      if (iablvisc == 2) pottempB=pottempA
   endif
   tau = 3.0*kinevisc + 0.5 ! static and constant tau also for ghost zones in case they are used

   call cpufinish(1)

#ifdef _CUDA
   if (ir == 0) call gpu_meminfo('time stepping')
#endif

! initialize pointers
   f1 => fA
   f2 => fB
   if (iablvisc == 2) then
      p1 => pottempA
      p2 => pottempB
   endif
   if (ntracer > 0) then
      t1 => tracerA
      t2 => tracerB
   endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (ir == 0) print *,'Start timestepping loop'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do it = nt0+1, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) then
         if (ir==0) write(*,'(a,i6,a,f10.2,a)')'Iteration: ',it,' Time: ',real(it)*p2l%time,' s.'
      endif

! start time step with f1,rho,u,v,w given

! External forcing
      external_forcing=0.0
      if (nturbines > 0)      call turbine_forcing(external_forcing, turbines, rho, u, v, w)
      if (iablvisc == 2)      call buoyancy_forcing(external_forcing,p1)

! [turbulence_df = turbulenceforcing(rho,u,v,w)]
      if (inflowturbulence)   call inflow_turbulence_forcing(rho,u,v,w,turbulence_ampl,it,nrturb)

! [f1,tau = post collision(f1,rho,u,v,w] (returns post collision density and tau for forcing)
      call postcoll(f1, tau, rho,u,v,w)


! [f1 = f1 + turbulence_df]
      if (inflowturbulence)   call inflow_turbulence_apply(f1,turbulence_df)

! [f1 = f1 + external forcing]
      if (nturbines > 0 .or. iablvisc==2) call forcings_apply(f1,external_forcing,rho,u,v,w)


! Bounce back boundary on fixed walls within the fluid
      if (lsolids) call solids(f1,lblanking)

! [f1 and f2 updated with boundary conditions]
      call boundarycond(f1,f2,uvel,t1,p1)

#ifdef MPI
      call mpi_halo_exchange_j(f1,nl)
#endif

! Drift of f1 returned in f2
      call drift(f2,f1)

! Swap such that f2,t2 becomes f1,t1 in next time step
      tmp   => f1
      f1  => f2
      f2 => tmp

! Compute updated macro variables and copy to halos for printing
#ifdef MPI
      call mpi_halo_exchange_j(f1,nl)
#endif
      call macrovars(rho,u,v,w,f1)

! pottemp advection returns updated pottemp in p2
      if (iablvisc == 2) then
#ifdef MPI
         call mpi_halo_exchange_j(p1,1)
#endif
         call heatflux(p2,p1)
         call advection(p2,p1,u,v,w,tau,1)
         pt_tmp => p1
         p1 => p2
         p2 => pt_tmp
      endif

! tracer advection returns updated tracer in t2
      if (ntracer > 0) then
#ifdef MPI
         call mpi_halo_exchange_j(t1,ntracer)
#endif
         call advection(t2,t1,u,v,w,tau,ntracer)
         tr_tmp => t1
         t1 => t2
         t2 => tr_tmp
      endif

! Diagnostics
#ifdef MPI
         call mpi_halo_exchange_j(t1,ntracer)
         call mpi_halo_exchange_j(p1,ntracer)
#endif
      call diag(itecout,it,rho,u,v,w,p1,t1,lblanking)

      call cpustart()
! Averaging for diagnostics similar to Asmuth paper for wind turbines
      if (laveturb .and. nturbines > 0) then
         if (avestart < it .and. it < avesave) call averaging_sec(u,v,w,.false.,turbines(1)%iradius)
         if (it == avesave)                    call averaging_sec(u,v,w,.true.,turbines(1)%iradius)
      endif

! Averaging for diagnostics
      if (laveraging) then
         if (avestart < it .and. it < avesave) call averaging_full(u,v,w,rho,p1,t1,lblanking,.false.)
         if (it == avesave)                    call averaging_full(u,v,w,rho,p1,t1,lblanking,.true.)
      endif

! Updating input turbulence matrix
      if (mod(it, nrturb) == 0 .and. it > 1 .and. inflowturbulence .and. ibnd==1) then
         call inflow_turbulence_update(uu,vv,ww,rr,nrturb,.false.)
      endif

! Save restart file
      if (mod(it,irestart) == 0)            call saverestart(it,f1,turbines(1:nturbines)%theta,uu,vv,ww,rr,p1,t1)
      call cpufinish(15)

      if (lmeasurements .and. mod(it,1000)==0) call predicted_measurements(u,v,w,it)

   enddo


   call cpustart()
   call saverestart(it-1,f1,turbines(1:nturbines)%theta,uu,vv,ww,rr,p1,t1)
   call cpufinish(16)
   if (ir==0) call cpuprint()

   if (allocated(uu       ))  deallocate(uu            )
   if (allocated(vv       ))  deallocate(vv            )
   if (allocated(ww       ))  deallocate(ww            )
   if (allocated(rr       ))  deallocate(rr            )
   if (allocated(tracerA  ))  deallocate(tracerA       )
   if (allocated(tracerB  ))  deallocate(tracerB       )

   call save_uvw(it,u,v,w)
   call testing(it,f1,f2)

#ifdef MPI
   call mpi_decomp_finalize()
#endif
end program LatticeBoltzmann

