module m_diag
contains
subroutine diag(filetype,it,rho,u,v,w,lblanking,Ti)
   use mod_dimensions
   use m_readinfile, only : iout, iprt1, iprt2, dprt,  nt1
   use m_vorticity
   use m_tecout
   use m_wtime
   implicit none
   integer, intent(in)   :: filetype
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
   character(len=10) cit

   character(len=200) :: variables='i,j,k,blanking,rho,u,v,w'
   real, allocatable    :: u_h(:,:,:)
   real, allocatable    :: v_h(:,:,:)
   real, allocatable    :: w_h(:,:,:)
   real, allocatable    :: rho_h(:,:,:)
   logical, allocatable :: lblanking_h(:,:,:)
   real, allocatable    :: Ti_h(:,:,:)
   integer, parameter :: icpu=14
   integer num_of_vars
   call cpustart()
   if ((mod(it, iout) == 0) .or. it == nt1 .or. ((it <= iprt1 .or. it >= iprt2).and. mod(it,dprt) == 0)) then

      if (filetype==0) then
         if (present(Ti)) then
            cit='F_AVERAGE'
            variables='i,j,k,blanking,rho,u,v,w,Ti'
            num_of_vars=9
         else
            write(cit,'(a1,i6.6)')'F',it
            variables='i,j,k,blanking,rho,u,v,w'
            num_of_vars=8
         endif

      elseif (filetype==1) then
         cit='_GRID'
         variables='i,j,k,blanking'
         num_of_vars=4

      elseif (filetype==2) then
         if (present(Ti)) then
            cit='_AVERAGE'
            variables='rho,u,v,w,Ti'
            num_of_vars=5
         else
            write(cit,'(i6.6)')it
            variables='rho,u,v,w'
            num_of_vars=4
         endif
      else
         print *,'diag filetype:',filetype
         stop 'invalid filetype in diag'
      endif

      allocate(u_h(nx,ny,nz))
      allocate(v_h(nx,ny,nz))
      allocate(w_h(nx,ny,nz))
      allocate(rho_h(nx,ny,nz))

      allocate(lblanking_h(0:nx+1,0:ny+1,0:nz+1))
      lblanking_h=lblanking

      u_h=u
      v_h=v
      w_h=w
      rho_h=rho

      if (present(Ti)) then
         allocate(Ti_h(nx,ny,nz))
         Ti_h=Ti
      endif

      if  (minval(rho_h) < 0.0) then
         print *,'iter=',it,'  minmaxrho=',minval(rho_h),' -- ',maxval(rho_h)
         print *,'iter=',it,'  minmaxloc=',minloc(rho_h),' -- ',maxloc(rho_h)
         stop 'Unstable simulation'
      endif

      if (present(Ti)) then
         call tecout(filetype,'tec'//trim(cit)//'.plt',it,trim(variables),num_of_vars,lblanking_h,rho_h,u_h,v_h,w_h,Ti_h)
      else
         call tecout(filetype,'tec'//trim(cit)//'.plt',it,trim(variables),num_of_vars,lblanking_h,rho_h,u_h,v_h,w_h)
      endif

      deallocate(u_h)
      deallocate(v_h)
      deallocate(w_h)
      deallocate(rho_h)
      if (allocated(lblanking_h)) deallocate(lblanking_h)
      if (allocated(Ti_h))        deallocate(Ti_h)
   endif

   call cpufinish(icpu)

end subroutine
end module
