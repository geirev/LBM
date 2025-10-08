module m_saverestart
contains
subroutine saverestart(it,f,theta,uu,vv,ww,rr)
   !use iso_fortran_env, only : int64
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile, only : inflowturbulence,nturbines,nrturb,lnodump
   implicit none
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
   integer iunit

   if (lnodump) return

   write(cit,'(i6.6)')it
   print '(a,a)',' saverestart:',cit
   if (inflowturbulence) then
      open(newunit=iunit,file='turbulence'//cit//'.uf',form="unformatted", status='unknown')
         uu_h=uu
         vv_h=vv
         ww_h=ww
         rr_h=rr
         write(iunit)ny,nz,nrturb,uu_h,vv_h,ww_h,rr_h
      close(iunit)
   endif

   if (nturbines > 0) then
      open(newunit=iunit,file='theta'//cit//'.dat')
         write(iunit,*)theta
      close(iunit)
   endif

   open(newunit=iunit,file='restart'//cit//'.uf',form="unformatted", status='replace')
      f_h=f
      write(iunit)nx,ny,nz,nl,f_h
   close(iunit)

end subroutine
end module

