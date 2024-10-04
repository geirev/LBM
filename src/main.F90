program LatticeBoltzmann
   use mod_dimensions
   use mod_D3Q27setup
   use m_readinfile
   use m_diag
   use m_averaging

   use m_airfoil
   use m_channel
   use m_cube
   use m_cylinder
   use m_disks
   use m_sphere
   use m_windfarm
   use m_channeldiag

   use m_uvelshear
   use m_boundarycond
   use m_bndbounceback
   use m_bndpressure
   use m_initurbulence

   use m_density
   use m_velocity
   use m_stress
   use m_turbineforcing
   use m_applyturbines
   use m_collisions
   use m_drift
   use m_fequil
   use m_readrestart
   use m_saverestart
   use m_wtime
   use, intrinsic :: omp_lib
   implicit none

! Main variables
   real    :: f(  0:nx+1,0:ny+1,0:nz+1,nl)   = 0.0  ! density function
   real    :: feq(0:nx+1,0:ny+1,0:nz+1,nl)   = 0.0  ! Maxwells equilibrium density function
   real    :: feqinflow(0:ny+1,0:nz+1,nl)    = 0.0  ! A scalar f used for inflow boundary conditions
   logical :: lblanking(nx,ny,nz)= .false.          ! blanking boundary

! Spatially dependent relaxation time
   real    :: tau(nx,ny,nz)      = 0.0              ! relaxation time scale

! Fluid variables
   real    :: u(nx,ny,nz)        = 0.0              ! x component of fluid velocity
   real    :: v(nx,ny,nz)        = 0.0              ! y component of fluid velocity
   real    :: w(nx,ny,nz)        = 0.0              ! z component of fluid velocity
   real    :: rho(nx,ny,nz)      = 0.0              ! fluid density
   real    :: sigma(3,3,nx,ny,nz)

! Stochastic input field on inflow boundary
   real uu(ny,nz,0:nrturb)
   real vv(ny,nz,0:nrturb)
   real ww(ny,nz,0:nrturb)
   real rr(ny,nz,0:nrturb)
   real inflowvar,inflowcor

! Turbine forcing
   real, allocatable  :: df(:,:,:,:,:)              ! Turbine forcing
   real uvel(nz)                                    ! vertical u-velocity profile

   integer :: it,i,j,l,k
   integer ip,jp,kp
   logical ex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading all input parameters
   call readinfile()
   if (nturbines > 0) allocate(df(-ieps:ieps,1:ny,1:nz,1:nl,nturbines))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define solid elements and walls
   select case(trim(experiment))
   case('cube')
      call cube(lblanking,nx/6,ny/2,nz/2,7)
   case('sphere')
      call sphere(lblanking,nx/6,ny/2,nz/2,10)
   case('cylinder')
      call cylinder(lblanking,nx/6,ny/2,25)
   case('windfarm')
      call windfarm(lblanking,kbnd)
   case('channel') ! Pouisille flow
      call channel(lblanking,nx/6,ny/2,25)
   case('airfoil')
      call airfoil(lblanking)
   case('disks')
      call disks(lblanking)
   case default
      print *,'invalid experiment',trim(experiment)
      stop
   end select

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initialization requires specification of u,v,w, and rho to compute feq
   call uvelshear(uvel)
! Initial turbulent flow
   inflowvar=0.001
   inflowcor=0.95
   tau=tauin

   if (nt0 == 0) then
      if (lpseudo) call initurbulence(uu,vv,ww,rr,rho,u,v,w,inflowcor,.true.)
      rho=rho0 + inflowvar*rho
      do k=1,nz
         u(:,:,k)=uvel(k)+inflowvar*u(:,:,k)
      enddo
      v=0.0+inflowvar*v
      w=0.0+inflowvar*w
      call diag(0,rho,u,v,w,lblanking)            ! Initial diagnostics

! Pressure gradient initialization for periodic boundaries with pressure drive.
!     if (ibnd==2) then
!        do i=1,nx
!           rho(i,:,:)=rho(i,:,:)+ rhoa - 2.0*rhoa*(real(i-1)/real(nx-1))
!        enddo
!     endif

! Final inititialization with equilibrium distribution from u,v,w, and rho
      call fequil(f,feq,rho,u,v,w,tau)            ! returns f for initialization

   else
! Restart from restart file
      call readrestart(nt0,f,theta,uu,vv,ww,rr)
      rho=density(f,lblanking)                    ! macro density
      u= velocity(f,rho,cxs,lblanking)            ! macro uvel
      v= velocity(f,rho,cys,lblanking)            ! macro vvel
      w= velocity(f,rho,czs,lblanking)            ! macro wvel
      call fequil(feq,f,rho,u,v,w,tau)            ! To get an initial value of tau for turbine forcing
   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation Main Loop
   open(99,file='checkforce.dat')
   do it = nt0+1, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) print '(a,i6,a,f10.2,a,a,f10.4)','Iteration:', it,' Time:',real(it)*p2l%time,' s'&
            ,' tau=',tau(nx/2,ny/2,nz/2)
                                                  ! start with f,rho,u,v,w
      call turbineforcing(df,rho,u,v,w,tau)       ! [u,v,w,df]         = turbineforcing(rho,u,v,w)
      call fequil(feq,f,rho,u,v,w,tau)            ! [feq,f=R(fneq),tau]= fequil(f,rho,u,v,w))
      call collisions(f,feq,tau)                  ! [feq=f]            = collisions(f,feq,tau)        f=f^eq + (1-1/tau) * R(f^neq)
      call applyturbines(feq,df,tau)              ! [feq=f]            = applyturbines(feq=f,df,tau)  f=f+df
      call boundarycond(feq,rho,u,v,w,rr,uu,vv,ww,it,inflowvar,uvel)  ! General boundary conditions
      call bndbounceback(feq,lblanking)           ! Bounce back boundary on fixed walls
      call drift(f,feq)                           ! Drift of feq returned in f
      rho=density(f,lblanking)                    ! macro density
      u= velocity(f,rho,cxs,lblanking)            ! macro uvel
      v= velocity(f,rho,cys,lblanking)            ! macro vvel
      w= velocity(f,rho,czs,lblanking)            ! macro wvel
      call diag(it,rho,u,v,w,lblanking)           ! Diagnostics

      ip=ipos(1)
      jp=jpos(1)
      kp=kpos(1)
      !write(99,'(i6,19f8.5,19f8.5)')it,u(ip:ip+36:2,jp-11,kp),w(ip:ip+36:2,jp-11,kp)
      write(99,'(i6,21f8.5,21f8.5)')it,u(ip-10:ip+10,jp-11,kp),w(ip-10:ip+10,jp-11,kp)



      if (avestart < it .and. it < avesave) call averaging(u,v,w,.false.,iradius)
      if (it == avesave)                    call averaging(u,v,w,.true.,iradius)
      if (mod(it, nrturb) == 0 .and. it > 1 .and. lpseudo) then
         print '(a,i6)','Recomputing inflow noise: it=',it
         call initurbulence(uu,vv,ww,rr,rho,u,v,w,inflowcor,.false.)
      endif
      if (mod(it,irestart) == 0)            call saverestart(it,f,theta,uu,vv,ww,rr)

   enddo
   close(99)

   call cpuprint()

   call saverestart(it-1,f,theta,uu,vv,ww,rr)

   select case(trim(experiment))
   case('channel')
      call channeldiag(f,rho,u,lblanking)
   end select

   if (allocated(df)) deallocate(df)

end program LatticeBoltzmann
