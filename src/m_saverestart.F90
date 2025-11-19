module m_saverestart
contains
subroutine saverestart(it,f,theta,uu,vv,ww,rr)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile, only : inflowturbulence,nturbines,nrturb,lnodump
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif
   implicit none
   integer, intent(in) :: it
   real,    intent(in) :: theta(nturbines)
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
#ifdef MPI
   character(len=4) ctile
#endif
   character(len=3)   ext
   character(len=5)   suffix
   character(len=10)  prefix
   character(len=10)  directory
   character(len=100) fname

   if (lnodump) return

! File names
#ifdef MPI
   write(ctile,'(i4.4)') mpi_rank
   suffix = '_' // trim(ctile)
#else
   suffix = ''
#endif
   ext='.uf'
   write(cit,'(i6.6)')it
   print '(4a)',' saverestart: tile=',trim(ctile),' iteration=',trim(cit)

   directory='restart/'
   call system('mkdir -p '//trim(directory))

   if (inflowturbulence) then
      prefix='turbulence'
      fname = trim(directory) // trim(prefix) // trim(suffix) // '_' // trim(cit) // trim(ext)
      open(newunit=iunit,file=trim(fname),form="unformatted", status='unknown')
         uu_h=uu
         vv_h=vv
         ww_h=ww
         rr_h=rr
         write(iunit)ny,nz,nrturb,uu_h,vv_h,ww_h,rr_h
      close(iunit)
   endif

   if (nturbines > 0) then
      prefix='theta'
      fname =  trim(directory) // trim(prefix) // trim(suffix) // '_' // trim(cit) // trim(ext)
      open(newunit=iunit,file=trim(fname),form="unformatted", status='replace')
         write(iunit)theta
      close(iunit)
   endif

   prefix='restart'
   fname =  trim(directory) // trim(prefix) // trim(suffix) // '_' // trim(cit) // trim(ext)
   open(newunit=iunit,file=trim(fname),form="unformatted", status='replace')
      f_h=f
      write(iunit)nx,ny,nz,nl,f_h
   close(iunit)

end subroutine
end module

