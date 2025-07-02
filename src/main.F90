program LatticeBoltzmann
   use mod_dimensions
   use mod_D3Q27setup
   use mod_shapiro
   use m_readinfile
   use m_assigncvel
   use m_diag
   use m_averaging
   use m_airfoil
   use m_city
   use m_cube
   use m_cylinder
   use m_sphere
   use m_channeldiag
   use m_uvelshear
   use m_boundarycond
   use m_bndpressure
   use m_inipert
   use m_initurbulence
   use m_macrovars
   use m_density
   use m_velocity
   use m_stress
   use m_solids
   use m_turbineforcing
   use m_turbulenceforcing
   use m_applyturbines
   use m_applyturbulence
   use m_collisions
   use m_drift
   use m_fequil
   use m_fequil3
   use m_fregularization
   use m_vreman
   use m_seedmanagement
   use m_readrestart
   use m_saverestart
   use m_wtime
   use m_set_random_seed2
   use, intrinsic :: omp_lib
   implicit none

   integer, parameter :: nshapiro=4
   real sh(0:nshapiro)


! Main variables
   real    :: f(nl,0:nx+1,0:ny+1,0:nz+1)     = 0.0  ! density function
   real    :: feq(nl,0:nx+1,0:ny+1,0:nz+1)   = 0.0  ! Maxwells equilibrium density function
   logical :: lblanking(0:nx+1,0:ny+1,0:nz+1)= .false.  ! blanking boundary and object grid points
   logical :: lsolids=.false.

! Spatially dependent relaxation time
   real    :: tau(nx,ny,nz)      = 0.0              ! relaxation time scale

! Fluid variables
   real    :: u(nx,ny,nz)        = 0.0              ! x component of fluid velocity
   real    :: v(nx,ny,nz)        = 0.0              ! y component of fluid velocity
   real    :: w(nx,ny,nz)        = 0.0              ! z component of fluid velocity
   real    :: rho(nx,ny,nz)      = 0.0              ! fluid density

! Stochastic input field on inflow boundary
   real uu(ny,nz,0:nrturb)
   real vv(ny,nz,0:nrturb)
   real ww(ny,nz,0:nrturb)
   real rr(ny,nz,0:nrturb)

! Turbine forcing
   real, allocatable  :: df(:,:,:,:,:)              ! Turbine forcing
   real, allocatable  :: turb_df(:,:,:,:)         ! Turbulence forcing
   real uvel(nz)                                    ! vertical u-velocity profile

   integer :: it
   integer ip,jp,kp

   logical, parameter :: runtest=.false.

   call set_random_seed2()

   !call assigncvel()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading all input parameters
   call readinfile()

   if (nturbines > 0) allocate(df(1:nl,-ieps:ieps,1:ny,1:nz,nturbines))
   if (lturb)         allocate(turb_df(1:nl,-ieps:ieps,1:ny,1:nz))
   !allocate(turb_df(1:nl,-ieps:ieps,1:ny,1:nz))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define solid elements and walls
   select case(trim(experiment))
   case('cube')
      call cube(lsolids,lblanking,nx/6,ny/2,nz/2,7)
   case('sphere')
      call sphere(lsolids,lblanking,nx/2,ny/2,nz/2,10)
   case('city')
      call city(lsolids,lblanking)
   case('cylinder')
      call cylinder(lsolids,lblanking,nx/2,ny/2,5)
   case('airfoil')
      call airfoil(lsolids,lblanking)
   end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization requires specification of u,v,w, and rho to compute feq
   call uvelshear(uvel)
   tau=tauin

! setting seed to seed.orig if file exist and nt0=0, otherwise generate new seed
   call seedmanagement(nt0)

! setting shapiro factors (no really used)
   call shfact(nshapiro,sh)

   if (nt0 == 0) then
! Intialization of macro variables
      call inipert(rho,u,v,w,uvel)

! Generate trubulence forcing fields
      if (lturb) call initurbulence(uu,vv,ww,rr,.true.)

! Outputs for testing
      if (runtest) then
         open(98,file='test.dat')
         ip=1; jp=jpos(1)-11; kp=kpos(1)
         ip=1; jp=ny/2; kp=nz/2
         write(98,'(a,100g13.5)')'000:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)
         !write(98,'(a,100g13.5)')'111:',uu(jp,kp,0:5),vv(jp,kp,0:5),ww(jp,kp,0:5),rr(jp,kp,0:5)
      endif

! Initial diagnostics
      call diag(0,rho,u,v,w,lblanking)


! Inititialization with equilibrium distribution from u,v,w, and rho
      call fequil3(feq,rho,u,v,w)
      call boundarycond(feq,rho,u,v,w,uvel)
      f=feq
      tau=tauin

   else
! Restart from restart file
      call readrestart(nt0,f,theta,uu,vv,ww,rr)
      call macrovars(rho,u,v,w,f,lblanking)

! To recover initial tau
      call fequil3(feq,rho,u,v,w)
      call fregularization(f, feq, u, v, w)
      call vreman(f,tau)
      f=feq+f
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation Main Loop
   do it = nt0+1, nt1
      if (runtest) write(98,'(a,100g13.5)')'001:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)
!      if ((mod(it, 10) == 0) .or. it == nt1) then
         write(*,'(a,i6,a,f10.2,a,3(a,f12.7))',advance='no')'Iteration:', it,                     &
                                                            ' Time:'  ,real(it)*p2l%time,' s,',   &
                                                            ' uinave:',sum(u(1,:,:))/real(ny*nz), &
                                                            ' rinave:',sum(rho(1,:,:))/real(ny*nz), &
                                                            ' tau:'   ,tau(nx/2,ny/2,nz/2)
!         do l=1,nl
!            write(*,'(a)',advance='no')'.'
!            call shfilt3D(nshapiro,sh,nshapiro,f(l,:,:,:),nx+2,ny+2,nz+2)
!         enddo
!         write(*,'(a)')'filtered'
         write(*,'(a)')'.'
!      endif

! start with f,rho,u,v,w
      if (runtest) then
         write(98,'(a,i7)')'it: ',it
         write(98,'(a,100g13.5)')'u02:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)
         write(98,'(a,100g13.5)')'f02:',f(1:10,ip,jp,kp)
      endif

! [u,v,w,df] = turbineforcing[rho,u,v,w]
      if (nturbines > 0) call turbineforcing(df,rho,u,v,w)

! [u,v,w,turb_df] = turbulenceforcing[rho,u,v,w]
      if (lturb) call turbulenceforcing(turb_df,rho,u,v,w,uu,vv,ww,it)

! [feq] = fequil3(rho,u,v,w] (returns equilibrium density)
      call fequil3(feq,rho,u,v,w)
      if (runtest) write(98,'(a,100g13.5)')'ini:',feq(1:10,ip,jp,kp)

! [f=Rneqf] = fregularization[f,feq,u,v,w] (input f is full f and returns reg. non-eq-density)
      call fregularization(f, feq, u, v, w)
      if (runtest) write(98,'(a,100g13.5)')'reg:',feq(1:10,ip,jp,kp)
      if (runtest) write(98,'(a,100g13.5)')'reg:',f(1:10,ip,jp,kp)

! [tau] = vreman[f] [f=Rneqf]
      call vreman(f,tau)
      if (runtest) write(98,'(a,100g13.5)')'vre:',tau(ip,jp,kp)

! [feq=f] = collisions(f,feq,tau)  f=f^eq + (1-1/tau) * R(f^neq)
      call collisions(f,feq,tau)
      if (runtest) write(98,'(a,100g13.5)')'col:',feq(1:10,ip,jp,kp)

! [feq=f] = applyturbines(feq,df,tau)  f=f+df
      if (nturbines > 0) call applyturbines(feq,df,tau)
      if (runtest) write(98,'(a,100g13.5)')'tur:',feq(1:10,ip,jp,kp)

! [feq=f] = applyturbulence(feq,turb_df,tau)  f=f+turb_df
      if (lturb) call applyturbulence(feq,turb_df,tau)
      if (runtest) write(98,'(a,100g13.5)')'tur:',feq(1:10,ip,jp,kp)

! Bounce back boundary on fixed walls within the fluid
      if (lsolids) call solids(feq,lblanking)

! General boundary conditions
      call boundarycond(feq,rho,u,v,w,uvel)
      if (runtest) write(98,'(a,100g13.5)')'bnd:',feq(1:10,ip,jp,kp)

! Drift of feq returned in f
      call drift(f,feq)
      if (runtest) write(98,'(a,100g13.5)')'dri:',f(1:10,ip,jp,kp)


! Compute updated macro variables
      call macrovars(rho,u,v,w,f,lblanking)
      write(98,'(a,100g13.5)')'u03:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)

! Diagnostics
      call diag(it,rho,u,v,w,lblanking)

! Averaging for diagnostics
      if (avestart < it .and. it < avesave) call averaging(u,v,w,.false.,iradius)
      if (it == avesave)                    call averaging(u,v,w,.true.,iradius)

! Updating input turbulence matrix
      if (mod(it, nrturb) == 0 .and. it > 1 .and. lturb .and. ibnd==1) then
         print '(a,i6)','Recomputing inflow noise: it=',it
         call initurbulence(uu,vv,ww,rr,.false.)
      endif

! Save restart file
      if (mod(it,irestart) == 0)            call saverestart(it,f,theta,uu,vv,ww,rr)
      if (runtest) write(98,'(a)')

   enddo

   if (runtest) close(98)

   call cpuprint()

   call saverestart(it-1,f,theta,uu,vv,ww,rr)

   select case(trim(experiment))
   case('channel')
      call channeldiag(f,rho,u,lblanking)
   end select

   if (allocated(df)) deallocate(df)

end program LatticeBoltzmann
