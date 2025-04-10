program LatticeBoltzmann
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
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
   use m_bndbounceback
   use m_closedbnd
   use m_bndpressure
   use m_inipert
   use m_initurbulence
   use m_macrovars
   use m_density
   use m_velocity
   use m_stress
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

! Main variables
   real    :: f(nl,0:nx+1,0:ny+1,0:nz+1)     = 0.0  ! density function
   real    :: feq(nl,0:nx+1,0:ny+1,0:nz+1)   = 0.0  ! Maxwells equilibrium density function
   logical :: lblanking(0:nx+1,0:ny+1,0:nz+1)= .false.  ! blanking boundary and object grid points

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

   logical, parameter :: forcecheck=.false.
   logical, parameter :: runtest=.false.

   call set_random_seed2()

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
      call cube(lblanking,nx/6,ny/2,nz/2,7)
   case('sphere')
      call sphere(lblanking,nx/2,ny/2,nz/2,10)
   case('city')
      call city(lblanking)
   case('cylinder')
      call cylinder(lblanking,nx/2,ny/2,5)
   case('airfoil')
      call airfoil(lblanking)
   end select

   call closedbnd(lblanking)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization requires specification of u,v,w, and rho to compute feq
   call uvelshear(uvel)
   tau=tauin

! setting seed to seed.orig if exist and nt0=0, otherwise generate new seed
   call seedmanagement(nt0)

   if (nt0 == 0) then
! Intialization of macro variables
      call inipert(rho,u,v,w,uvel)

! Generate trubulence forcing fields
      if (lturb) call initurbulence(uu,vv,ww,rr,.true.)

!      if (ibnd==2) then ! Pressure gradient initialization for periodic boundaries with pressure drive.
!         stop 'ibnd=2 does not work I think'
!         do i=1,nx
!            rho(i,:,:)=rho0 + inflowstd*rho(i,:,:) +  rhoa - 2.0*rhoa*(real(i-1)/real(nx-1))
!         enddo
!         u=0.0; v=0.0; w=0.0
!      endif

! Outputs for testing
      if (runtest) then
         open(98,file='test.dat')
         ip=2; jp=jpos(1)-11; kp=kpos(1)
         write(98,'(a,100e19.10)')'000:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)
         write(98,'(a,100e19.10)')'111:',uu(jp,kp,0:5),vv(jp,kp,0:5),ww(jp,kp,0:5),rr(jp,kp,0:5)
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

! to recover initial tau
      call fequil3(feq,rho,u,v,w)
      call fregularization(f, feq, rho, u, v, w)
      call vreman(f,tau)
      f=feq+f
   endif

   if (forcecheck) then
      open(99,file='checkforce.dat')
      ip=ipos(1); jp=jpos(1); kp=kpos(1)
      write(99,'(i6,21f9.5,21f9.5)')it,u(ip-10:ip+10,jp-11,kp),w(ip-10:ip+10,jp-11,kp)
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation Main Loop
   do it = nt0+1, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) print '(a,i6,a,f10.2,a,a,f10.4)','Iteration:', it,' Time:',real(it)*p2l%time,' s,',&
                                           ' uinave:',sum(u(1,:,:))/real(ny*nz)

! start with f,rho,u,v,w
      if (runtest) then
         ip=2; jp=jpos(1)-11; kp=kpos(1)
         write(98,'(a,i7)')'ite:',it
         write(98,'(a,100e19.10)')'ini:',u(ip,jp,kp),v(ip,jp,kp),w(ip,jp,kp),rho(ip,jp,kp)
         write(98,'(a,100e19.10)')'ini:',f(1:25,ip,jp,kp)
      endif

! [u,v,w,df] = turbineforcing[rho,u,v,w]
      if (nturbines > 0) call turbineforcing(df,rho,u,v,w)

! [u,v,w,turb_df] = turbineforcing[rho,u,v,w]
      if (lturb) call turbulenceforcing(turb_df,rho,u,v,w,uu,vv,ww,it)

! [feq] = fequil3(rho,u,v,w] (returns equilibrium density)
      call fequil3(feq,rho,u,v,w)
      if (runtest) write(98,'(a,100e19.10)')'ini:',feq(1:25,ip,jp,kp)

! [f=Rneqf] = fregularization[f,feq,u,v,w] (input f is full f and returns reg. non-eq-density)
      call fregularization(f, feq, rho, u, v, w)
      if (runtest) write(98,'(a,100e19.10)')'reg:',feq(1:25,ip,jp,kp)
      if (runtest) write(98,'(a,100e19.10)')'reg:',f(1:25,ip,jp,kp)

! [tau] = vreman[f] [f=Rneqf]
      call vreman(f,tau)
      if (runtest) write(98,'(a,100e19.10)')'vre:',tau(ip,jp,kp)

! [feq=f] = collisions(f,feq,tau)  f=f^eq + (1-1/tau) * R(f^neq)
      call collisions(f,feq,tau)
      if (runtest) write(98,'(a,100e19.10)')'col:',feq(1:25,ip,jp,kp)

! [feq=f] = applyturbines(feq,df,tau)  f=f+df
      if (nturbines > 0) call applyturbines(feq,df,tau)
      if (runtest) write(98,'(a,100e19.10)')'tur:',feq(1:25,ip,jp,kp)

! [feq=f] = applyturbulence(feq,turb_df,tau)  f=f+turb_df
      if (lturb) call applyturbulence(feq,turb_df,tau)
      if (runtest) write(98,'(a,100e19.10)')'tur:',feq(1:25,ip,jp,kp)

! General boundary conditions
      call boundarycond(feq,rho,u,v,w,uvel)
      if (runtest) write(98,'(a,100e19.10)')'bnd:',feq(1:25,ip,jp,kp)

! Bounce back boundary on fixed walls
      call bndbounceback(feq,lblanking)
      if (runtest) write(98,'(a,100e19.10)')'bon:',feq(1:25,ip,jp,kp)

! Drift of feq returned in f
      call drift(f,feq)
      if (runtest) write(98,'(a,100e19.10)')'dri:',f(1:25,ip,jp,kp)
      if (runtest) write(98,'(a)')

! Compute updated macro variables
      call macrovars(rho,u,v,w,f,lblanking)

! Diagnostics
      call diag(it,rho,u,v,w,lblanking)

      if (forcecheck) then
         ip=ipos(1); jp=jpos(1); kp=kpos(1)
         write(99,'(i6,21f9.5,21f9.5)')it,u(ip-10:ip+10,jp-11,kp),w(ip-10:ip+10,jp-11,kp)
      endif

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

   enddo

   if (forcecheck) close(99)
   if (runtest) close(98)

   call cpuprint()

   call saverestart(it-1,f,theta,uu,vv,ww,rr)

   select case(trim(experiment))
   case('channel')
      call channeldiag(f,rho,u,lblanking)
   end select

   if (allocated(df)) deallocate(df)

end program LatticeBoltzmann
