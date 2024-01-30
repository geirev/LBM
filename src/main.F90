program LatticeBoltzmann
   use mod_dimensions
   use mod_D3Q27setup
   use m_initialization
   use m_airfoil
   use m_bndapply
   use m_channel
   use m_collisions
   use m_cube
   use m_cylinder
   use m_defbnd
   use m_density
   use m_disks
   use m_drift
   use m_fequil
   use m_gnuplot
   use m_gnuplotuf
   use m_gnuplotuv
   use m_readinfile
   use m_readrestart
   use m_saverestart
   use m_sphere
   use m_tecplot3D
   use m_velocity
   use m_vorticity
   use m_wtime
   use, intrinsic :: omp_lib
   implicit none


! Main variables
   real    :: f(nx,ny,nz,nl)     = 0.0        ! density function
   real    :: feq(nx,ny,nz,nl)   = 0.0        ! Maxwells equilibrium density function
   real    :: fbnd(nx,ny,nz,nl)  = 0.0        ! Boundary conditions
   logical :: blanking(nx,ny,nz) = .false.    ! blanking boundary

! Fluid variables
   real    :: u(nx,ny,nz)        = 0.0        ! x component of fluid velocity
   real    :: v(nx,ny,nz)        = 0.0        ! y component of fluid velocity
   real    :: w(nx,ny,nz)        = 0.0        ! z component of fluid velocity
   real    :: rho(nx,ny,nz)      = 0.0        ! fluid density
   real    :: speed(nx,ny,nz)    = 0.0        ! absolute velocity
   real    :: vortx(nx,ny,nz)    = 0.0        ! fluid vorticity x-component
   real    :: vorty(nx,ny,nz)    = 0.0        ! fluid vorticity y-component
   real    :: vortz(nx,ny,nz)    = 0.0        ! fluid vorticity z-component
   real    :: vort(nx,ny,nz)     = 0.0        ! absolute value of vorticity

   integer :: i,j,k,l,it,ia,ja,ib,jb
   character(len=6) cit

   call readinfile

   select case(trim(experiment))
   case('cube')
      call cube(blanking,nx/6,ny/2,nz/2,7)
      f=initialization(rho0,3.0,0.5,0.0,blanking)
   case('sphere')
      call sphere(blanking,nx/6,ny/2,nz/2,10)
      f=initialization(rho0,12.0,1.5,0.0,blanking)
   case('cylinder')
      call cylinder(blanking,nx/6,ny/2,25)
      f=initialization(rho0,1.5,0.5,0.0,blanking)
   case('channel')
      call channel(blanking,nx/6,ny/2,25)
      f=initialization(rho0,1.5,0.5,0.0,blanking)
   case('airfoil')
      call airfoil(blanking)
      f=initialization(rho0,2.5,0.0,0.0,blanking)
   case('disks')
      call disks(blanking)
      f=initialization(rho0,2.5,0.0,0.0,blanking)
   case default
      print *,'invalid experiment',trim(experiment)
      stop
   end select

! Restart from restart file
   if (nt0 > 1) then
      write(cit,'(i6.6)')nt0
      call readrestart('restart'//cit//'.uf',f)
   endif

! Calculate fluid variables and plotting of initial conditions
   if (nt0 == 0) then
      write(cit,'(i6.6)')nt0
      rho=density(f,blanking)
      u= velocity(f,rho,cxs,blanking)
      v= velocity(f,rho,cys,blanking)
      w= velocity(f,rho,czs,blanking)
      speed = sqrt(u*u + v*v + w*w)
      call  vorticity(u,v,w,vortx,vorty,vortz,vort,blanking)

!     call gnuplotuf('vortz'//cit//'.uf',vortz(:,:,nz/2),nx,ny)
!     call gnuplotuv('uv'//cit//'.dat',u(:,:,nz/2),v(:,:,nz/2),nx,ny)
      call tecplot3D('tec'//cit//'.dat',rho,u,v,w,speed,vortx,vorty,vortz,vort,blanking)
   endif

! Simulation Main Loop
   fbnd=0.0
   do it = nt0, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) print '(a,i6)','iteration:', it

! Drift
      cpu0=wtime(); call cpu_time(start)
      call drift(f,feq)
      cpu1=wtime(); call cpu_time(finish)
      cputime(1)=cputime(1)+finish-start
      waltime(1)=waltime(1)+cpu1-cpu0

! Set reflective boundaries
      cpu0=wtime(); call cpu_time(start)
      call defbnd(fbnd,f,blanking)
      cpu1=wtime(); call cpu_time(finish)
      cputime(2)=cputime(2)+finish-start
      waltime(2)=waltime(2)+cpu1-cpu0

! Calculate equlibrium distribution
      cpu0=wtime(); call cpu_time(start)
      rho=density(f,blanking)
      cpu1=wtime(); call cpu_time(finish)
      cputime(3)=cputime(3)+finish-start
      waltime(3)=waltime(3)+cpu1-cpu0

      cpu0=wtime(); call cpu_time(start)
      u= velocity(f,rho,cxs,blanking)
      v= velocity(f,rho,cys,blanking)
      w= velocity(f,rho,czs,blanking)
      cpu1=wtime(); call cpu_time(finish)
      cputime(4)=cputime(4)+finish-start
      waltime(4)=waltime(4)+cpu1-cpu0

      cpu0=wtime(); call cpu_time(start)
      call fequil(feq,rho,u,v,w)
      cpu1=wtime(); call cpu_time(finish)
      cputime(5)=cputime(5)+finish-start
      waltime(5)=waltime(5)+cpu1-cpu0

! Apply collisions
      cpu0=wtime(); call cpu_time(start)
      call collisions(f,feq,tau)
      cpu1=wtime(); call cpu_time(finish)
      cputime(6)=cputime(6)+finish-start
      waltime(6)=waltime(6)+cpu1-cpu0

! Apply boundary
      cpu0=wtime(); call cpu_time(start)
      call bndapply(f,fbnd,blanking)
      cpu1=wtime(); call cpu_time(finish)
      cputime(7)=cputime(7)+finish-start
      waltime(7)=waltime(7)+cpu1-cpu0

! Plots dv/dx-du/dy
      cpu0=wtime(); call cpu_time(start)
      if ((mod(it, iout) == 0) .or. it == nt1) then
         print '(a)','Dumping diagnostics...'
         speed = sqrt(u*u + v*v + w*w)
         call vorticity(u,v,w,vortx,vorty,vortz,vort,blanking)
         write(cit,'(i6.6)')it
!         call gnuplot('vort'//cit//'.dat',vorty,nx,ny)
!         call gnuplotuf('vortz'//cit//'.uf',vortz,nx,ny)
!         call gnuplotuv('uv'//cit//'.dat',u(:,:,nz/2),v(:,:,nz/2),nx,ny)
         call tecplot3D('tec'//cit//'.dat',rho,u,v,w,speed,vortx,vorty,vortz,vort,blanking)
      endif
      cpu1=wtime(); call cpu_time(finish)
      cputime(8)=cputime(8)+finish-start
      waltime(8)=waltime(8)+cpu1-cpu0

   enddo
   print '(tr22,3a)','cputime  ','walltime  ','speedup   '
   do l=1,8
      print '(tr10,a9,3f10.4)',cpuname(l),cputime(l),waltime(l),cputime(l)/waltime(l)
   enddo
   print '(tr10,a9,3f10.4)','summary  ',sum(cputime(1:8)),sum(waltime(1:8)),sum(cputime(1:8))/sum(waltime(1:8))
! Dump restart file
   write(cit,'(i6.6)')it-1
   call saverestart('restart'//cit//'.uf',f)

end program LatticeBoltzmann
