program LatticeBoltzmann
   use mod_dimensions
   use mod_D3Q27setup
   use m_set_random_seed2
   use m_pseudo2D
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

   use m_ibndinflow
   use m_boundarycond
   use m_bndbounceback
   use m_bndpressure

   use m_density
   use m_velocity
   use m_stress
   use m_turbineforcing
   use m_applyturbines
   use m_collisions
   use m_drift
   use m_BGKequil
   use m_HRRequil
   use m_feqscalar
   use m_fhrrscalar
   use m_readrestart
   use m_saverestart
   use m_wtime
   use, intrinsic :: omp_lib
   implicit none



! Main variables
   real    :: f(  0:nx+1,0:ny+1,0:nz+1,nl)   = 0.0  ! density function
   real    :: feq(0:nx+1,0:ny+1,0:nz+1,nl)   = 0.0  ! Maxwells equilibrium density function
   real    :: feqscal(0:nz+1,nl)             = 0.0  ! A scalar f used for inflow boundary conditions
   logical :: lblanking(nx,ny,nz)= .false.          ! blanking boundary

! Spatially dependent relaxation time
   real    :: tau(nx,ny,nz)      = 0.0          ! x component of fluid velocity
! Fluid variables
   real    :: u(nx,ny,nz)        = 0.0          ! x component of fluid velocity
   real    :: v(nx,ny,nz)        = 0.0          ! y component of fluid velocity
   real    :: w(nx,ny,nz)        = 0.0          ! z component of fluid velocity
   real    :: rho(nx,ny,nz)      = 0.0          ! fluid density
   real    :: sigma(3,3,nx,ny,nz)

! Turbine forcing
   real, allocatable  :: df(:,:,:,:,:)            ! Turbine forcing

   integer :: it,i,j,l,k
   real cor1,cor2,dx,dy,dir
   real xlength,ylength
   integer n1,n2
   real width,mu,grad,cs2,dw
   logical ex


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reading all input parameters
   call readinfile()
   if (nturbines > 0) allocate(df(-ieps:ieps,1:ny,1:nz,1:nl,nturbines))

   tau(:,:,:)=tauin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Sample 2D pseudo random fields
   call set_random_seed2
   cor1=20.0/sqrt(3.0)
   cor2=20.0/sqrt(3.0)
   dir=0.0
   xlength=real(nx-1)
   ylength=real(ny-1)
   dx=1.0
   dy=1.0
   n1=nx
   n2=ny

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
   u=uini
   v=0.0
   w=0.0
   rho=rho0

! Adding smooth density perturbations
   if (lpseudo) then
     call pseudo2D(rho,nx,ny,nz,cor1,cor2,dx,dy,n1,n2,dir,.true.)
     call pseudo2D(u,nx,ny,nz,cor1,cor2,dx,dy,n1,n2,dir,.true.)
     call pseudo2D(v,nx,ny,nz,cor1,cor2,dx,dy,n1,n2,dir,.true.)
     call pseudo2D(w,nx,ny,nz,cor1,cor2,dx,dy,n1,n2,dir,.true.)
     do k=2,nz
        rho(:,:,k)=0.1*rho(:,:,k-1)+sqrt(1.0-0.1**2)*rho(:,:,k)
        u(:,:,k)=0.1*u(:,:,k-1)+sqrt(1.0-0.1**2)*u(:,:,k)
        v(:,:,k)=0.1*v(:,:,k-1)+sqrt(1.0-0.1**2)*v(:,:,k)
        w(:,:,k)=0.1*w(:,:,k-1)+sqrt(1.0-0.1**2)*w(:,:,k)
     enddo
     rho=rho0 + 0.0001*rho
     u=uini+0.0001*u
     v=0.0+0.0001*v
     w=0.0+0.0001*w
   endif

   if (ibnd==1) then
! Constant inflow boundary velocity read from file and reinitialization of u(nx,ny,nz)
! Specification of constant-in-time inflow boundary distribution
      call ibndinflow(feqscal,u)

   elseif (ibnd==2) then
! Pressure gradient initialization for periodic boundaries with pressure drive.
      do i=1,nx
         rho(i,:,:)=rho(i,:,:)+ rhoa - 2.0*rhoa*(real(i-1)/real(nx-1))
      enddo
   endif

! Final inititialization with equilibrium distribution from u,v,w, and rho
   if (coll=='HRR') then
      call HRRequil(f,feq,rho,u,v,w,tau)  ! returns f for initialization
   elseif (coll=='BGK') then
      call BGKequil(f,feq,rho,u,v,w)
   endif

! Restart from restart file
   if (nt0 > 1) call readrestart(f,nt0,theta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Simulation Main Loop
   do it = nt0, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) print '(a,i6,a,f10.2,a)','Iteration:', it,' Time:',real(it)*p2l%time,' s'

      rho=density(f,lblanking)                    ! macro density
      u= velocity(f,rho,cxs,lblanking)            ! macro uvel
      v= velocity(f,rho,cys,lblanking)            ! macro vvel
      w= velocity(f,rho,czs,lblanking)            ! macro wvel
      call diag(it,rho,u,v,w,lblanking)           ! Diagnostics
      if (coll=='HRR') then
         call HRRequil(feq,f,rho,u,v,w,tau)       ! f is input, returns feq, R(fneq) in f, and tau
      elseif (coll=='BGK') then
         call BGKequil(feq,f,rho,u,v,w)           ! Calculate equlibrium distribution
      endif
      call turbineforcing(df,feq,rho,u,v,w)       ! define forcing df from each turbine
      call collisions(f,feq,tau)                  ! Apply collisions, returns feq
      call applyturbines(feq,df)                  ! operates on f stored in feq
      call boundarycond(feq,rho,u,v,w,feqscal)    ! General boundary conditions
      call bndbounceback(feq,lblanking)           ! Bounce back boundary on fixed walls
      call drift(f,feq)                           ! Drift of feq returned in f
      if (avestart < it .and. it < avesave) call averaging(u,v,w,.false.,iradius)
      if (it == avesave)                    call averaging(u,v,w,.true.,iradius)
      if (mod(it,irestart) == 0)            call saverestart(f,it,theta)

   enddo

   call cpuprint()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call saverestart(f,it-1,theta)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   select case(trim(experiment))
   case('channel')
      cs2=1.0/3.0
      rho=density(f,lblanking)
      u= velocity(f,rho,cxs,lblanking)
      width=real(ny-2)
      mu=cs2*(tauin-0.5)
      grad=cs2*2.0*rhoa/real(nx)
      dw=width/real(ny-1)
      open(10,file='channeluvel.dat')
      do j=1,ny
         write(10,'(i4,2f13.5)')j,u(nx/2,j,2),&
                    (1.0/(2.0*rho0*mu))*grad*((width/2.0)**2-(dw*(real(j)-1.0)-width/2.0)**2)
      enddo
      close(10)
   end select

   if (allocated(df)) deallocate(df)

end program LatticeBoltzmann
