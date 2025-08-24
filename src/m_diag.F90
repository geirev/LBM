module m_diag
contains
subroutine diag(it,rho,u,v,w,lblanking,Ti)
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
   real, optional,  intent(in)   :: Ti(nx,ny,nz)
#ifdef _CUDA
   attributes(device) :: rho
   attributes(device) :: u
   attributes(device) :: v
   attributes(device) :: w
   attributes(device) :: lblanking
   attributes(device) :: Ti
#endif
   character(len=6) cit

   integer, parameter :: num_of_fvars=31
   character(len=200) :: variables='i,j,k,blanking,rho,u,v,w'
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
   real, allocatable    :: Ti_h(:,:,:)
   integer, parameter :: icpu=14
   integer num_of_vars
   real tmp
   integer i,j,k
   call cpustart()
   if ((mod(it, iout) == 0) .or. it == nt1 .or. ((it <= iprt1 .or. it >= iprt2).and. mod(it,dprt) == 0)) then
      if (present(Ti)) then
         cit='averag'
         variables='i,j,k,blanking,rho,u,v,w,Ti'
      else
         write(cit,'(i6.6)')it
         variables='i,j,k,blanking,rho,u,v,w'
      endif

      allocate(u_h(nx,ny,nz))
      allocate(v_h(nx,ny,nz))
      allocate(w_h(nx,ny,nz))
      allocate(rho_h(nx,ny,nz))
      allocate(lblanking_h(0:nx+1,0:ny+1,0:nz+1))

      u_h=u
      v_h=v
      w_h=w
      rho_h=rho
      lblanking_h=lblanking
      num_of_vars=8

      if (present(Ti)) then
         allocate(Ti_h(nx,ny,nz))
         num_of_vars=num_of_vars+1
         Ti_h=Ti
      endif

      if (minval(rho_h) < 0.0) then
         print *,'iter=',it,'  minmaxrho=',minval(rho_h),' -- ',maxval(rho_h)
         print *,'iter=',it,'  minmaxloc=',minloc(rho_h),' -- ',maxloc(rho_h)
         stop 'Unstable simulation'
      endif
      call tecout('tec'//cit//'.plt',it,trim(variables),num_of_vars,lblanking_h,rho_h,u_h,v_h,w_h,Ti_h)
      deallocate(u_h)
      deallocate(v_h)
      deallocate(w_h)
      deallocate(rho_h)
      deallocate(lblanking_h)
   endif

   call cpufinish(icpu)

end subroutine
end module
