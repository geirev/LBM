program LatticeBoltzmann
   use mod_dimensions
   use mod_D3Q27setup
   use m_set_random_seed2
   use m_pseudo2D
   use m_readinfile
   use m_diag

   use m_airfoil
   use m_channel
   use m_cube
   use m_cylinder
   use m_disks
   use m_sphere
   use m_windfarm

   use m_boundarycond
   use m_bndbounceback
   use m_bndpressure
   use m_phys2lattice

   use m_density
   use m_velocity
   use m_collisions
   use m_drift
   use m_fequil
   use m_feqscalar
   use m_readrestart
   use m_saverestart
   use m_wtime
   use, intrinsic :: omp_lib
   implicit none



! Main variables
   real    :: f(  0:nx+1,0:ny+1,0:nz+1,nl)   = 0.0  ! density function
   real    :: feq(0:nx+1,0:ny+1,0:nz+1,nl)   = 0.0  ! Maxwells equilibrium density function
   real    :: feqscal(nl)                    = 0.0  ! A scalar f used for boundary conditions
   logical :: lblanking(nx,ny,nz)= .false.          ! blanking boundary

! Fluid variables
   real    :: u(nx,ny,nz)        = 0.0          ! x component of fluid velocity
   real    :: v(nx,ny,nz)        = 0.0          ! y component of fluid velocity
   real    :: w(nx,ny,nz)        = 0.0          ! z component of fluid velocity
   real    :: rho(nx,ny,nz)      = 0.0          ! fluid density

   integer :: it,i,j
   character(len=6) cit

   real cor1,cor2,dx,dy,dir
   real xlength,ylength
   integer n1,n2


   real width,mu,grad,cs2,dw

   call readinfile

   call phys2lattice
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
      call windfarm(lblanking,ny/2,nz/2,5)
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

! Initialization
   u=uini
   v=0.0
   w=0.0
   rho=rho0
   if (lpseudo) then
     call pseudo2D(rho,nx,ny,nz,cor1,cor2,dx,dy,n1,n2,dir,.true.)
     rho=rho0 + 0.1*rho
   endif

   if (ibnd==1) then
! Constant inflow boundary distribution.
      call feqscalar(feqscal,rho0,u(1,1,1),v(1,1,1),w(1,1,1))

   elseif (ibnd==2) then
! Pressure gradient initialization for periodic boundaries with pressure drive.
      do i=1,nx
         rho(i,:,:)=rho(i,:,:)+ rhoa - 2.0*rhoa*(real(i-1)/real(nx-1))
      enddo
      call feqscalar(feqscal,rho0,u(1,1,1),v(1,1,1),w(1,1,1))
   endif

! Inititialization with equilibrium distribution
   call fequil(f,rho,u,v,w)

! Restart from restart file
   if (nt0 > 1) then
      write(cit,'(i6.6)')nt0
      call readrestart('restart'//cit//'.uf',f)
   endif

! Simulation Main Loop
   do it = nt0, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) print '(a,i6)','iteration:', it

      rho=density(f,lblanking)                    ! macro density
      u= velocity(f,rho,cxs,lblanking)            ! macro uvel
      v= velocity(f,rho,cys,lblanking)            ! macro vvel
      w= velocity(f,rho,czs,lblanking)            ! macro wvel

      call fequil(feq,rho,u,v,w)                  ! Calculate equlibrium distribution

      call diag(it,rho,u,v,w,lblanking)           ! Diagnostics

      call collisions(f,feq,tau)                  ! Apply collisions, returns feq

      call boundarycond(feq,rho,u,v,w,feqscal)    ! General boundary conditions

      call bndbounceback(feq,lblanking)           ! Bounce back boundary on fixed walls

      call drift(f,feq)                           ! Drift

   enddo

   call cpuprint()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Dump restart file
   write(cit,'(i6.6)')it-1
   call saverestart('restart'//cit//'.uf',f)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   select case(trim(experiment))
   case('channel')
      cs2=1.0/3.0
      rho=density(f,lblanking)
      u= velocity(f,rho,cxs,lblanking)
      width=real(ny-2)
      mu=cs2*(tau-0.5)
      grad=cs2*2.0*rhoa/real(nx)
      dw=width/real(ny-1)
      open(10,file='uvel.dat')
      do j=1,ny
         write(10,'(i4,2f13.5)')j,u(nx/2,j,2),&
                    (1.0/(2.0*rho0*mu))*grad*((width/2.0)**2-(dw*(real(j)-1.0)-width/2.0)**2)
      enddo
      close(10)
   end select

end program LatticeBoltzmann
