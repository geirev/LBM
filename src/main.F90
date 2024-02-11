program LatticeBoltzmann
   use mod_dimensions
   use mod_D3Q27setup
   use m_pseudo2D
   use m_set_random_seed2
   use m_tecout
   use m_pdfout
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
   use m_feqscalar
   use m_gnuplot
   use m_gnuplotuf
   use m_gnuplotuv
   use m_readinfile
   use m_readrestart
   use m_saverestart
   use m_sphere
!   use m_tecplot3D
   use m_velocity
   use m_vorticity
   use m_wtime
   use, intrinsic :: omp_lib
   implicit none


! Main variables
   real    :: f(nx,ny,nz,nl)     = 0.0        ! density function
   real    :: feq(nx,ny,nz,nl)   = 0.0        ! Maxwells equilibrium density function
   real    :: fbnd(nx,ny,nz,nl)  = 0.0        ! Boundary conditions
   logical :: lblanking(nx,ny,nz)= .false.    ! blanking boundary
   real    :: feqscal(nl)        = 0.0        ! A scalar f used for boundary conditions

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

   character(len=200)::tecplot_fvars='i,j,k,blanking,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27'
   integer, parameter :: num_of_fvars=31

   character(len=200) :: tecplot_variables='i,j,k,blanking,rho,u,v,w,vel,vortx,vorty,vortz,vort'
   integer, parameter :: num_of_variables=13


   real cor1,cor2,dx,dy,dir
   real xlength,ylength
   integer n1,n2

   call readinfile

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Set a variable random seed
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

! Sample 2D pseudo random fields
   call pseudo2D(rho,nx,ny,nz,cor1,cor2,dx,dy,n1,n2,dir,.true.)


   select case(trim(experiment))
   case('cube')
      call cube(lblanking,nx/6,ny/2,nz/2,7)
   case('sphere')
      call sphere(lblanking,nx/6,ny/2,nz/2,10)
   case('cylinder')
      call cylinder(lblanking,nx/6,ny/2,25)
   case('channel')
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
!   call initialization(f,lblanking) ! old initialization
   u=0.12
   v=0.0
   w=0.0
!   call random_number(rho)
   rho=rho0 + 0.1*rho
   call fequil(f,rho,u,v,w)
   call feqscalar(feqscal,rho0,u(1,1,1),v(1,1,1),w(1,1,1))

! Restart from restart file
   if (nt0 > 1) then
      write(cit,'(i6.6)')nt0
      call readrestart('restart'//cit//'.uf',f)
   endif

! Simulation Main Loop
   fbnd=0.0
   do it = nt0, nt1
      if ((mod(it, 10) == 0) .or. it == nt1) print '(a,i6)','iteration:', it

! Inflow boundary
      do l=1,nl
      do k=1,nz
      do j=1,ny
         f(1,j,k,l)=feqscal(l)
         f(nx,j,k,l)=f(nx-1,j,k,l)
      enddo
      enddo
      enddo

! macro variables
      cpu0=wtime(); call cpu_time(start)
      rho=density(f,lblanking)
      cpu1=wtime(); call cpu_time(finish)
      cputime(3)=cputime(3)+finish-start
      waltime(3)=waltime(3)+cpu1-cpu0

      cpu0=wtime(); call cpu_time(start)
      u= velocity(f,rho,cxs,lblanking)
      v= velocity(f,rho,cys,lblanking)
      w= velocity(f,rho,czs,lblanking)
!      print '(a,3g13.5)','dens : ',rho(200,50,2)
!      print '(a,3g13.5)','vel u: ',u(200,50,2)
!      print '(a,3g13.5)','vel v: ',v(200,50,2)
!      print '(a,3g13.5)','vel w: ',w(200,50,2)
      cpu1=wtime(); call cpu_time(finish)
      cputime(4)=cputime(4)+finish-start
      waltime(4)=waltime(4)+cpu1-cpu0

! Calculate equlibrium distribution
      cpu0=wtime(); call cpu_time(start)
      call fequil(feq,rho,u,v,w)
!      print '(a,27g13.5)','feq  : ',feq(200,50,2,1:27)
      cpu1=wtime(); call cpu_time(finish)
      cputime(5)=cputime(5)+finish-start
      waltime(5)=waltime(5)+cpu1-cpu0

! Diagnostics
      cpu0=wtime(); call cpu_time(start)
      if ((mod(it, iout) == 0) .or. it == nt1) then
         print '(a)','Dumping diagnostics...'
         speed = sqrt(u*u + v*v + w*w)
         call vorticity(u,v,w,vortx,vorty,vortz,vort,lblanking)
         write(cit,'(i6.6)')it
         call tecout('tec'//cit//'.plt',trim(tecplot_variables),num_of_variables,lblanking,rho,u,v,w,speed,vortx,vorty,vortz,vort)
      endif
!      if ((mod(it, ifout) == 0)) then
!         write(cit,'(i6.6)')it
!         call pdfout('pdf'//cit//'.plt',trim(tecplot_fvars),num_of_fvars,lblanking,f)
!      endif
      cpu1=wtime(); call cpu_time(finish)
      cputime(8)=cputime(8)+finish-start
      waltime(8)=waltime(8)+cpu1-cpu0

! Apply collisions
      cpu0=wtime(); call cpu_time(start)
      call collisions(f,feq,tau)
      cpu1=wtime(); call cpu_time(finish)
      cputime(6)=cputime(6)+finish-start
      waltime(6)=waltime(6)+cpu1-cpu0

! Drift
      cpu0=wtime(); call cpu_time(start)
      call drift(f,feq)
      cpu1=wtime(); call cpu_time(finish)
      cputime(1)=cputime(1)+finish-start
      waltime(1)=waltime(1)+cpu1-cpu0

! Set reflective boundaries
      cpu0=wtime(); call cpu_time(start)
      call defbnd(fbnd,f,lblanking)
      cpu1=wtime(); call cpu_time(finish)
      cputime(2)=cputime(2)+finish-start
      waltime(2)=waltime(2)+cpu1-cpu0

! Apply boundary
      cpu0=wtime(); call cpu_time(start)
      call bndapply(f,fbnd,lblanking)
      cpu1=wtime(); call cpu_time(finish)
      cputime(7)=cputime(7)+finish-start
      waltime(7)=waltime(7)+cpu1-cpu0


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
