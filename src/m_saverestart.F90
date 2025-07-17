module m_saverestart
contains
subroutine saverestart(it,f,theta,uu,vv,ww,rr)
   use mod_dimensions
   use m_readinfile, only : inflowturbulence,nturbines
   integer, intent(in) :: it
   real,    intent(in) :: theta
   real,    intent(in) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in) :: uu(ny,nz,0:nrturb)
   real,    intent(in) :: vv(ny,nz,0:nrturb)
   real,    intent(in) :: ww(ny,nz,0:nrturb)
   real,    intent(in) :: rr(ny,nz,0:nrturb)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
   attributes(device) :: rr
#endif
   real :: f_h(nl,0:nx+1,0:ny+1,0:nz+1)
   real :: uu_h(ny,nz,0:nrturb)
   real :: vv_h(ny,nz,0:nrturb)
   real :: ww_h(ny,nz,0:nrturb)
   real :: rr_h(ny,nz,0:nrturb)

   character(len=6) cit
   integer :: irec

   write(cit,'(i6.6)')it
   print '(a,a)',' saverestart:',cit
   if (inflowturbulence) then
      inquire(iolength=irec)ny,nz,nrturb,uu_h,vv_h,ww_h,rr_h
      open(10,file='turbulence'//cit//'.uf',form="unformatted", access="direct", recl=irec)
         uu_h=uu
         vv_h=vv
         ww_h=ww
         rr_h=rr
         write(10,rec=1)ny,nz,nrturb,uu_h,vv_h,ww_h,rr_h
      close(10)
   endif

   if (nturbines > 0) then
      open(10,file='theta'//cit//'.dat')
         write(10,*)theta
      close(10)
   endif

   inquire(iolength=irec) nx,ny,nz,nl,f_h
   open(10,file='restart'//cit//'.uf',form="unformatted", access="direct", recl=irec)
      f_h=f
      write(10,rec=1)nx,ny,nz,nl,f_h
   close(10)

end subroutine
end module

