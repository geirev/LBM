module m_diag
contains
subroutine diag(it,rho,u,v,w,lblanking)
   use mod_dimensions
   use m_readinfile
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
   logical, intent(in)   :: lblanking(nx,ny,nz)
   character(len=6) cit

   character(len=200)::tecplot_fvars='i,j,k,blanking,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27'
   integer, parameter :: num_of_fvars=31
   character(len=200) :: tecplot_variables='i,j,k,blanking,rho,u,v,w,vel,vortx,vorty,vortz,vort'
   integer, parameter :: num_of_variables=13
   real    :: speed(nx,ny,nz)    = 0.0        ! absolute velocity
   real    :: vortx(nx,ny,nz)    = 0.0        ! fluid vorticity x-component
   real    :: vorty(nx,ny,nz)    = 0.0        ! fluid vorticity y-component
   real    :: vortz(nx,ny,nz)    = 0.0        ! fluid vorticity z-component
   real    :: vort(nx,ny,nz)     = 0.0        ! absolute value of vorticity
   integer, parameter :: icpu=8





   call cpustart()
   if ((mod(it, iout) == 0) .or. it == nt1) then
      if (minval(rho) < 0.0) then
         print *,'iter=',it,'  minmaxrho=',minval(rho),' -- ',maxval(rho)
         stop 'Unstable simulation'
      endif

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
   call cpufinish(icpu)

end subroutine
end module
