program LatticeBoltzmann
   use mod_dimensions
   use m_initialization
   use m_cylinder
   use m_airfoil
   use m_density
   use m_velocity
   use m_drift
   use m_gnuplot
   use m_gnuplotuf
   use m_gnuplotuv
   use m_defbnd
   use m_fequil
   use m_collisions
   use m_wtime
   use, intrinsic :: omp_lib
   implicit none

! Simulation parameters
   integer, parameter :: nt   = 8000        ! 8000       ! number of timesteps
   real,    parameter :: rho0 = 100.0      ! average density
   real,    parameter :: tau  = 0.62       ! collision timescale 0.6
   logical, parameter :: dbg  =.false.     ! print diagnostics if tru

! Main variables
   real    :: f(nx,ny,nl)     = 0.0        ! density function
   real    :: feq(nx,ny,nl)   = 0.0        ! Maxwells equilibrium density function
   real    :: fbnd(nx,ny,nl)  = 0.0        ! Boundary conditions
   logical :: blanking(nx,ny) = .false.    ! blanking boundary

! Fluid variables
   real    :: u(nx,ny)        = 0.0        ! x component of fluid velocity
   real    :: v(nx,ny)        = 0.0        ! y component of fluid velocity
   real    :: rho(nx,ny)      = 0.0        ! fluid density
   real    :: vorticity(nx,ny)= 0.0        ! fluid vorticity

! Lattice speeds D2Q9 and weights
   integer :: cxs(1:nl) = [0, 0, 1, 1, 1, 0, -1, -1, -1]
   integer :: cys(1:nl) = [0, 1, 1, 0, -1, -1, -1, 0, 1]
   real    :: weights(1:nl) = [4.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0, 1.0/9.0, 1.0/36.0]

   integer :: i,j,l,it,ia,ja,ib,jb
   character(len=4) cit
   character(len=20) :: experiment='airfoil'
!   character(len=20) :: experiment='cylinder'



   select case(trim(experiment))
   case('cylinder')
      f=initialization(rho0,1.5,0.5)
      call cylinder(blanking,nx/6,ny/2,25)
   case('airfoil')
      f=initialization(rho0,2.5,0.0)
      call airfoil(blanking)
   case default
      print *,'invalid experiment',trim(experiment)
      stop
   end select


! Calculate fluid variables and plotting of initial conditions
   rho=density(f)
   u= velocity(f,rho,cxs)
   v= velocity(f,rho,cys)
   vorticity = (cshift(v,1,1) - cshift(v,-1,1)) - (cshift(u,1,2) - cshift(u,-1,2))
   where (blanking)
      u=0.0
      v=0.0
      vorticity=0.0
   endwhere
   it=0
   write(cit,'(i4.4)')it
!   call gnuplot('vort0000.dat',u,nx,ny)
   call gnuplotuf('vort0000.uf',vorticity,nx,ny)
   call gnuplotuv('uv'//cit//'.dat',u,v,nx,ny)

! Simulation Main Loop
   fbnd=0.0
   do it = 1, nt
      if ((mod(it, 1000) == 0) .or. it == Nt) print '(a,i4)','iteration:', it

! Drift
      cpu0=wtime(); call cpu_time(start)
      call drift(f,feq,cxs,cys)
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
      rho=density(f)
      cpu1=wtime(); call cpu_time(finish)
      cputime(3)=cputime(3)+finish-start
      waltime(3)=waltime(3)+cpu1-cpu0

      cpu0=wtime(); call cpu_time(start)
      u= velocity(f,rho,cxs)
      v= velocity(f,rho,cys)
      cpu1=wtime(); call cpu_time(finish)
      cputime(4)=cputime(4)+finish-start
      waltime(4)=waltime(4)+cpu1-cpu0

      cpu0=wtime(); call cpu_time(start)
      call fequil(feq,rho,u,v,weights,cxs,cys)
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
      do j=1,ny
      do i=1,nx
         if (blanking(i,j)) f(i,j,1:nl)=fbnd(i,j,1:nl)
      enddo
      enddo
      cpu1=wtime(); call cpu_time(finish)
      cputime(7)=cputime(7)+finish-start
      waltime(7)=waltime(7)+cpu1-cpu0

! Plots dv/dx-du/dy
      cpu0=wtime(); call cpu_time(start)
      if ((mod(it, 10) == 0) .or. it == Nt) then
!        vorticity = sqrt(u*u + v*v)
!!$OMP PARALLEL DO PRIVATE(i,j,ia,ib,ja,jb) SHARED(u,v)
!          do j=1,ny
!             ja=mod(j-2,ny)+1
!             jb=mod(j,ny)+1
!             do i=1,nx
!                ia=mod(i-2,nx)+1
!                ib=mod(i,nx)+1
!                vorticity(i,j) = (v(ib,j) - v(ia,j))   -   (u(i,jb) - u(i,ja))
!             enddo
!          enddo
!!$OMP END PARALLEL DO

         vorticity = (cshift(v,1,1) - cshift(v,-1,1)) - (cshift(u,1,2) - cshift(u,-1,2))
         where (blanking)
            u=0.0
            v=0.0
            vorticity=0.0
         endwhere
         write(cit,'(i4.4)')it
!         call gnuplot('vort'//cit//'.dat',vorticity,nx,ny)
         call gnuplotuf('vort'//cit//'.uf',vorticity,nx,ny)
         call gnuplotuv('uv'//cit//'.dat',u,v,nx,ny)
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


end program LatticeBoltzmann
