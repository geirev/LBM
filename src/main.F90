program LatticeBoltzmann
#ifdef _CUDA
   use cudafor
#endif
   use mod_dimensions
   use mod_D3Q27setup
   use mod_shapiro
   use m_readinfile
   use m_hermite_polynomials
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
   use m_initurbulence
   use m_macrovars
   use m_density
   use m_velocity
   use m_solids
   use m_turbineforcing
   use m_turbulenceforcing
   use m_applyturbines
   use m_applyturbulence
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
   use m_set_random_seed2
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

! Turbine forcing
   real, allocatable  :: df(:,:,:,:,:)              ! Turbine forcing
   real, allocatable  :: turb_df(:,:,:,:)         ! Turbulence forcing
#ifdef _CUDA
   attributes(device) :: df
   attributes(device) :: turb_df
#endif

   real :: elevation(nx,ny)=0.0
   integer i,j,k,l
   integer :: it
   integer ip,jp,kp
   integer :: istat

   logical, parameter :: runtest=.true.

#ifdef _CUDA
   print*, "Running GPU version!"
#else
   print*, "Running CPU version!"
#endif

   call set_random_seed2()
   call hermite_polynomials()

   !call assigncvel()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading all input parameters
   call readinfile()

   if (nturbines > 0) allocate(df(1:nl,-ieps:ieps,1:ny,1:nz,nturbines))
   if (lturb)         allocate(turb_df(1:nl,-ieps:ieps,1:ny,1:nz))
   !allocate(turb_df(1:nl,-ieps:ieps,1:ny,1:nz))

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
   case('cube')
!      call cube(lsolids,lblanking,nx/6,ny/2,nz/2,7)
   case('sphere')
!      call sphere(lsolids,lblanking,nx/2,ny/2,nz/2,10)
   case('city')
!      call city(lsolids,lblanking)
   case('cylinder')
      call cylinder(lsolids,lblanking)
      lblanking_h=lblanking
      do j=1,ny
      do i=1,nx
         if (lblanking_h(i,j,1)) then
            elevation(i,j)=nz
         endif
      enddo
      enddo
      call tecfld('elevation',nx,ny,1,elevation)
   case('airfoil')
!      call airfoil(lsolids,lblanking)
   end select

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
      if (lturb) call initurbulence(uu,vv,ww,rr,.true.)

! Outputs for testing
!      if (runtest) then
!         open(98,file='test.dat')
!         ip=1; jp=jpos(1)-11; kp=kpos(1)
!         ip=1; jp=ny/2; kp=nz/2
!         write(98,'(a,100g15.7)')'000:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)
!      endif

! Initial diagnostics
      call diag(0,rho,u,v,w,lblanking)


! Inititialization with equilibrium distribution from u,v,w, and rho
      call fequil3(feq,rho,u,v,w, A2, A3, vel)
!      write(98,'(a,100g15.7)')'fXX:',feq(1:10,ip,jp,kp)
      call boundarycond(feq,rho,u,v,w,uvel_d)
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k = 0, nz+1
      do j = 0, ny+1
      do i = 0, nx+1
         f(:,i,j,k) = feq(:,i,j,k)
      end do
      end do
      end do
!      write(98,'(a,100g15.7)')'f00:',f(1:10,ip,jp,kp)

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k = 1, nz
      do j = 1, ny
      do i = 1, nx
         tau(i,j,k) = tauin
      end do
      end do
      end do

   else
! Restart from restart file
      call readrestart(nt0,f,theta,uu,vv,ww,rr)
      call macrovars(rho,u,v,w,f,lblanking)

! To recover initial tau
      call fequil3(feq,rho,u,v,w, A2, A3, vel)
      call regularization(f, feq, u, v, w, A2, A3, vel)
      call vreman(f, tau, eddyvisc ,Bbeta ,alphamag ,alpha ,beta)

#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k = 0, nz+1
      do j = 0, ny+1
      do i = 0, nx+1
         f(:,i,j,k)=feq(:,i,j,k)+f(:,i,j,k)
      end do
      end do
      end do

   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation Main Loop
   do it = nt0+1, nt1
      if (it>=25000) iout=2
!      if (runtest) write(98,'(a,100g15.7)')'001:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)
      if ((mod(it, 10) == 0) .or. it == nt1) then
!          istat=cudaMemcpy(tmptau, tau(nx/2,ny/2,nz/2), 8, cudaMemcpyDeviceToHost)
!         print '(a)','x'
!          tmptau = work_h(nx/2,ny/2,nz/2)
!         print '(a)','Y'
!         work_h=u;   tmpu   = sum(work_h(1,:,:))/real(ny*nz)
!         print '(a)','Z'
!         work_h=rho; tmprho = sum(work_h(1,:,:))/real(ny*nz)
!         print '(a)','A'
         write(*,'(a,i6,a,f10.2,a,3(a,f12.7))',advance='no')'Iteration:', it,                     &
                                                            ' Time:'  ,real(it)*p2l%time,' s,'
!                                                            ' uinave:',tmpu,                      &
!                                                            ' rinave:',tmprho,                    &
                                                             !' tau:'   ,tmptau
          !deallocate(work_h)
!         do l=1,nl
!            write(*,'(a)',advance='no')'.'
!            call shfilt3D(nshapiro,sh,nshapiro,f(l,:,:,:),nx+2,ny+2,nz+2)
!         enddo
!         write(*,'(a)')'filtered'
         write(*,'(a)')'.'
      endif

! start with f,rho,u,v,w
!      if (runtest) then
!         write(98,'(a,i7)')'it: ',it
!         write(98,'(a,100g15.7)')'u02:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)
!         write(98,'(a,100g15.7)')'f02:',f(1:10,ip,jp,kp)
!      endif

! [u,v,w,df] = turbineforcing[rho,u,v,w]
      if (nturbines > 0) call turbineforcing(df,rho,u,v,w)

! [u,v,w,turb_df] = turbulenceforcing[rho,u,v,w]
      if (lturb) call turbulenceforcing(turb_df,rho,u,v,w,uu,vv,ww,it)

! [feq] = fequil3(rho,u,v,w] (returns equilibrium density)
      call fequil3(feq,rho,u,v,w, A2, A3, vel)
#ifdef _CUDA
     ! istat = cudaDeviceSynchronize()
#endif
!      if (runtest) write(98,'(a,100g15.7)')'ini:',feq(1:10,ip,jp,kp)

! [f=Rneqf] = regularization[f,feq,u,v,w] (input f is full f and returns reg. non-eq-density)
      call regularization(f, feq, u, v, w, A2, A3, vel)
#ifdef _CUDA
     ! istat = cudaDeviceSynchronize()
#endif

!      if (runtest) write(98,'(a,100g15.7)')'reg:',feq(1:10,ip,jp,kp)
!      if (runtest) write(98,'(a,100g15.7)')'reg:',f(1:10,ip,jp,kp)

! [tau] = vreman[f] [f=Rneqf]
      call vreman(f, tau, eddyvisc ,Bbeta ,alphamag ,alpha ,beta)
!      if (runtest) write(98,'(a,100g15.7)')'vre:',tau(ip,jp,kp)

! [feq=f] = collisions(f,feq,tau)  f=f^eq + (1-1/tau) * R(f^neq)
      call collisions(f,feq,tau)
!      if (runtest) write(98,'(a,100g15.7)')'col:',feq(1:10,ip,jp,kp)

! [feq=f] = applyturbines(feq,df,tau)  f=f+df
      if (nturbines > 0) call applyturbines(feq,df,tau)
!      if (runtest) write(98,'(a,100g15.7)')'tur:',feq(1:10,ip,jp,kp)

! [feq=f] = applyturbulence(feq,turb_df,tau)  f=f+turb_df
      if (lturb) call applyturbulence(feq,turb_df,tau)
!      if (runtest) write(98,'(a,100g15.7)')'tur:',feq(1:10,ip,jp,kp)

! Bounce back boundary on fixed walls within the fluid
      if (lsolids) call solids(feq,lblanking)

! General boundary conditions
      call boundarycond(feq,rho,u,v,w,uvel_d)
!      if (runtest) write(98,'(a,100g15.7)')'bnd:',feq(1:10,ip,jp,kp)

! Drift of feq returned in f
      call drift(f,feq)
!      if (runtest) write(98,'(a,100g15.7)')'dri:',f(1:10,ip,jp,kp)


! Compute updated macro variables
      call macrovars(rho,u,v,w,f,lblanking)
!      if (runtest) write(98,'(a,100g15.7)')'u03:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)

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

   if (allocated(df)) deallocate(df)

end program LatticeBoltzmann
