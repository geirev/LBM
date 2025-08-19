module m_diag
contains
subroutine diag(it,rho,u,v,w,lblanking)
   use mod_dimensions
   use m_readinfile, only : iout, iprt1, iprt2, dprt, lprtmin,  nt1
   use m_vorticity
   use m_tecout
   use m_pdfout
   use m_wtime
   implicit none
   integer, intent(in)   :: it
   real,    intent(in)   :: rho(nx,ny,nz)
   real,    intent(in)   :: u(nx,ny,nz)
   real,    intent(in)   :: v(nx,ny,nz)
   real,    intent(in)   :: w(nx,ny,nz)
   logical, intent(in)   :: lblanking(0:nx+1,0:ny+1,0:nz+1)
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: lblanking
#endif
   character(len=6) cit

!   character(len=200)::tecplot_fvars='i,j,k,blanking,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27'
   integer, parameter :: num_of_fvars=31
!   character(len=200) :: tecplot_maxvar='i,j,k,x,y,z,blanking,rho,u,v,w,vel,vortx,vorty,vortz,vort'
   character(len=200) :: tecplot_minvar='i,j,k,blanking,rho,u,v,w'
   real    :: speed(nx,ny,nz)           ! absolute velocity
   real    :: vortx(nx,ny,nz)           ! fluid vorticity x-component
   real    :: vorty(nx,ny,nz)           ! fluid vorticity y-component
   real    :: vortz(nx,ny,nz)           ! fluid vorticity z-component
   real    :: vort(nx,ny,nz)            ! absolute value of vorticity
#ifdef _CUDA
   attributes(device) :: speed
   attributes(device) :: vortx
   attributes(device) :: vorty
   attributes(device) :: vortz
   attributes(device) :: vort
#endif
   real, allocatable    :: u_h(:,:,:)
   real, allocatable    :: v_h(:,:,:)
   real, allocatable    :: w_h(:,:,:)
   real, allocatable    :: rho_h(:,:,:)
   logical, allocatable :: lblanking_h(:,:,:)
   integer, parameter :: icpu=14
   integer num_of_vars
   real tmp
   integer i,j,k
   call cpustart()
   if ((mod(it, iout) == 0) .or. it == nt1 .or. ((it <= iprt1 .or. it >= iprt2).and. mod(it,dprt) == 0)) then



      if (.not.lprtmin) then
         call vorticity(u,v,w,vortx,vorty,vortz,vort,lblanking)
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#else
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(speed,u,v,w)
#endif
         do k=1,nz
         do j=1,ny
         do i=1,nx
            tmp=u(i,j,k)*u(i,j,k) + v(i,j,k)*v(i,j,k) + w(i,j,k)*w(i,j,k)
            speed(i,j,k)=sqrt(tmp)
         enddo
         enddo
         enddo
#ifndef _CUDA
!$OMP END PARALLEL DO
#endif
      endif

      write(cit,'(i6.6)')it
      if (.not.lprtmin) then
        !num_of_vars=16
        !call tecout('tec'//cit//'.plt',it,trim(tecplot_maxvar),num_of_vars,lblanking,rho,u,v,w,&
        !            speed,vortx,vorty,vortz,vort)
      else
         allocate(u_h(nx,ny,nz))
         allocate(v_h(nx,ny,nz))
         allocate(w_h(nx,ny,nz))
         allocate(rho_h(nx,ny,nz))
         allocate(lblanking_h(0:nx+1,0:ny+1,0:nz+1))
         num_of_vars=8
         u_h=u
         v_h=v
         w_h=w
         rho_h=rho
         lblanking_h=lblanking

         if (minval(rho_h) < 0.0) then
            print *,'iter=',it,'  minmaxrho=',minval(rho_h),' -- ',maxval(rho_h)
            print *,'iter=',it,'  minmaxloc=',minloc(rho_h),' -- ',maxloc(rho_h)
            stop 'Unstable simulation'
         endif
         call tecout('tec'//cit//'.plt',it,trim(tecplot_minvar),num_of_vars,lblanking_h,rho_h,u_h,v_h,w_h)
         deallocate(u_h)
         deallocate(v_h)
         deallocate(w_h)
         deallocate(rho_h)
         deallocate(lblanking_h)
      endif
   endif

   call cpufinish(icpu)

end subroutine
end module
