program LatticeBoltzmann
#ifdef _CUDA
   use cudafor
#endif
   use m_rhotest
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
   use m_uvelshear
   use m_boundarycond
   use m_inipert
   use m_macrovars
   use m_density
   use m_velocity
   use m_solids
   use m_turbineforcing
   use m_applyturbulence
   use m_initurbulence
   use m_turbulenceforcing
   use m_applyturbines
   use m_collisions
   use m_drift
   use m_fequil3
   use m_regularization
   use m_vreman
   use m_seedmanagement
   use m_readrestart
   use m_saverestart
   use m_tecfld
   use m_wtime
   use m_set_random_seed3
   use, intrinsic :: omp_lib
   implicit none

   integer, parameter :: nshapiro=4
   real sh(0:nshapiro)


! Main variables
   real    :: f(nl,0:nx+1,0:ny+1,0:nz+1)            ! density function
   real    :: feq(nl,0:nx+1,0:ny+1,0:nz+1)          ! Maxwells equilibrium density function
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: feq
#endif

   logical :: lblanking(0:nx+1,0:ny+1,0:nz+1)       ! blanking boundary and object grid points
#ifdef _CUDA
   attributes(device) :: lblanking
#endif
   logical, dimension(:,:,:), allocatable :: lblanking_h ! Host blanking variable

   real uvel(nz)
   real, dimension(:), allocatable :: uvel_d        ! vertical u-velocity profile host version
#ifdef _CUDA
   attributes(device) :: uvel_d
#endif

   logical :: lsolids=.false.

! Spatially dependent relaxation time
   real    :: tau(nx,ny,nz)                         ! relaxation time scale
#ifdef _CUDA
   attributes(device) :: tau
#endif
   real, allocatable :: work_h(:,:,:)
   real tmptau,tmpu,tmprho

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

! Hermite coefficients fequil and regularization
   real :: vel(3,nx,ny,nz)
   real :: A2(3,3,nx,ny,nz)
   real :: A3(3,3,3,nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: vel
   attributes(device) :: A2
   attributes(device) :: A3
#endif

! Vreman arrays
   real eddyvisc(nx,ny,nz)
   real Bbeta(nx,ny,nz)
   real alphamag(nx,ny,nz)
   real alpha(3,3,nx,ny,nz)
   real beta(3,3,nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: alpha
   attributes(device) :: beta
   attributes(device) :: eddyvisc
   attributes(device) :: Bbeta
   attributes(device) :: alphamag
#endif

! Stochastic input field on inflow boundary
   real uu(ny,nz,0:nrturb)
   real vv(ny,nz,0:nrturb)
   real ww(ny,nz,0:nrturb)
   real rr(ny,nz,0:nrturb)
#ifdef _CUDA
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
   attributes(device) :: rr
#endif


   real, allocatable :: myturb(:,:,:)
   real :: elevation(nx,ny)=0.0
   integer i,j,k,l
   integer :: it
   integer ip,jp,kp

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call cpustart()
   call hermite_polynomials()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading all input parameters
   call readinfile()

   if (nturbines > 0)      call init_turbines()
   if (inflowturbulence)   call init_turbulenceforcing()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
   allocate(lblanking_h(0:nx,0:ny,0:nz))
   select case(trim(experiment))
   case('city')
      call city(lsolids,lblanking)
   case('cylinder')
      call cylinder(lsolids,lblanking)
   case('airfoil')
!      call airfoil(lsolids,lblanking)
       stop 'needs fix airfoil routine for gpu'
   end select

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
      if (inflowturbulence) call initurbulence(uu,vv,ww,rr,.true.)

! Initial diagnostics
      call diag(0,rho,u,v,w,lblanking)


! Inititialization with equilibrium distribution from u,v,w, and rho
      call fequil3(feq,rho,u,v,w, A2, A3, vel,it, nt1)
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
      do k=1,nz
      do j=1,ny
      do i=1,nx
         tau(i,j,k) = tauin
      enddo
      enddo
      enddo

   else
! Restart from restart file
      call readrestart(nt0,f,theta,uu,vv,ww,rr)
      call macrovars(rho,u,v,w,f)

! To recover initial tau
      call fequil3(feq,rho,u,v,w, A2, A3, vel,it,nt1)
      call regularization(f, feq, u, v, w, A2, A3, vel,it,nt1)
      call vreman(f, tau, eddyvisc ,Bbeta ,alphamag ,alpha ,beta, it,nt1)

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k=0,nz+1
      do j=0,ny+1
      do i=0,nx+1
         do l=1,nl
            f(:,i,j,k)=feq(:,i,j,k)+f(:,i,j,k)
         enddo
      enddo
      enddo
      enddo

   endif
   call cpufinish(1)

   print *,'Start timestepping loop'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation Main Loop
   do it = nt0+1, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) then
!          istat=cudaMemcpy(tmptau, tau(nx/2,ny/2,nz/2), 8, cudaMemcpyDeviceToHost)
!         print '(a)','x'
!          tmptau = work_h(nx/2,ny/2,nz/2)
!         print '(a)','Y'
!         work_h=u;   tmpu   = sum(work_h(1,:,:))/real(ny*nz)
!         print '(a)','Z'
!         work_h=rho; tmprho = sum(work_h(1,:,:))/real(ny*nz)
!         print '(a)','A'
         write(*,'(a,i6,a,f10.2,a,3(a,f12.7))',advance='yes')'Iteration:', it,                     &
                                                            ' Time:'  ,real(it)*p2l%time,' s.'
!                                                            ' uinave:',tmpu,                      &
!                                                            ' rinave:',tmprho,                    &
                                                             !' tau:'   ,tmptau
          !deallocate(work_h)
      endif

! start with f,rho,u,v,w


! [u,v,w,turbine_df] = turbineforcing[rho,u,v,w]
      if (nturbines > 0)      call turbineforcing(rho,u,v,w,it,nt1)

! [u,v,w,turbulence_df] = turbulenceforcing[rho,u,v,w]
      if (inflowturbulence)   call turbulenceforcing(rho,u,v,w,uu,vv,ww,turbulence_ampl,it,nt1)

! [feq] = fequil3(rho,u,v,w] (returns equilibrium density)
      call fequil3(feq,rho,u,v,w, A2, A3, vel,it, nt1);               if (debug) call rhotest(feq,rho,'fequil')

! [f=Rneqf] = regularization[f,feq,u,v,w] (input f is full f and returns reg. non-eq-density)
      call regularization(f, feq, u, v, w, A2, A3, vel,it, nt1)


! [tau] = vreman[f] [f=Rneqf]
      call vreman(f, tau, eddyvisc ,Bbeta ,alphamag ,alpha ,beta, it, nt1)

! [feq=f] = collisions(f,feq,tau)  f=f^eq + (1-1/tau) * R(f^neq)
      call collisions(f,feq,tau);                                     if (debug) call rhotest(feq,rho,'collisions')

! [feq=f] = applyturbines(feq,turbine_df,tau)  f=f+turbine_df
      if (nturbines > 0)      call applyturbines(feq,turbine_df,tau)

! [feq=f] = applyturbulence(feq,turbulence_df,tau)  f=f+turb_df
      if (inflowturbulence)   call applyturbulence(feq,turbulence_df,tau);  if (debug) call rhotest(feq,rho,'applyturbulence')

! Bounce back boundary on fixed walls within the fluid
      if (lsolids) call solids(feq,lblanking);                        if (debug) call rhotest(feq,rho,'solids')

! General boundary conditions
      call boundarycond(feq,uvel_d);                                  if (debug) call rhotest(feq,rho,'boundarycond')

! Drift of feq returned in f
      call drift(f,feq);                                              if (debug) call rhotest(f,rho,'drift')


! Compute updated macro variables
      call macrovars(rho,u,v,w,f);                                    if (debug) call rhotest(f,rho,'macrovars')

! Diagnostics
      call diag(it,rho,u,v,w,lblanking);                              if (debug) call rhotest(f,rho,'diag')

      call cpustart()
! Averaging for diagnostics
      if (avestart < it .and. it < avesave) call averaging(u,v,w,.false.,iradius)
      if (it == avesave)                    call averaging(u,v,w,.true.,iradius)

! Updating input turbulence matrix
      if (mod(it, nrturb) == 0 .and. it > 1 .and. inflowturbulence .and. ibnd==1) then
         print '(a,i6)','Recomputing inflow noise: it=',it
         call initurbulence(uu,vv,ww,rr,.false.)
      endif

! Save restart file
      if (mod(it,irestart) == 0)            call saverestart(it,f,theta,uu,vv,ww,rr)
      call cpufinish(15)

   enddo


   call cpustart()
   call saverestart(it-1,f,theta,uu,vv,ww,rr)
   call cpufinish(16)
   call cpuprint()

end program LatticeBoltzmann
