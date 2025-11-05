module m_readrestart
contains
subroutine readrestart(it,f,theta,uu,vv,ww,rr)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile, only : inflowturbulence,nturbines,nrturb
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif
   implicit none
   integer, intent(in)  :: it
   real,    intent(out) :: theta
   real,    intent(out) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(out) :: uu(ny,nz,0:nrturb)
   real,    intent(out) :: vv(ny,nz,0:nrturb)
   real,    intent(out) :: ww(ny,nz,0:nrturb)
   real,    intent(out) :: rr(ny,nz,0:nrturb)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
   attributes(device) :: rr
#endif
   real, allocatable  :: f_h(:,:,:,:)
   real, allocatable  :: uu_h(:,:,:)
   real, allocatable  :: vv_h(:,:,:)
   real, allocatable  :: ww_h(:,:,:)
   real, allocatable  :: rr_h(:,:,:)

   logical ex
   integer :: i,j,k,l,n
   character(len=6) cit
   integer iunit
#ifdef MPI
   character(len=4) ctile
#endif
   character(len=3) ext
   character(len=5) suffix
   character(len=10) prefix
   character(len=10)  directory
   character(len=100) fname

   allocate(f_h(nl,0:nx+1,0:ny+1,0:nz+1))
   allocate(uu_h(ny,nz,0:nrturb)        )
   allocate(vv_h(ny,nz,0:nrturb)        )
   allocate(ww_h(ny,nz,0:nrturb)        )
   allocate(rr_h(ny,nz,0:nrturb)        )

! File names
#ifdef MPI
   write(ctile,'(i4.4)') mpi_rank
   suffix = '_' // trim(ctile)
#else
   suffix = ''
#endif
   ext='.uf'
   write(cit,'(i6.6)')it

   directory='restart/'
   call system('mkdir -p '//trim(directory))

   if (inflowturbulence) then
      prefix='turbulence'
      fname =  trim(directory) // trim(prefix) // trim(suffix) // '_' // trim(cit) // trim(ext)
      inquire(file=trim(fname),exist=ex)
      print '(3a)','reading: ',trim(fname)
      if (ex) then
         open(newunit=iunit,file=trim(fname),form="unformatted", status='old')
            read(iunit,err=998)j,k,l
            rewind(iunit)
            if ((j==ny).and.(k==nz).and.(l==nrturb)) then
               read(iunit,err=998)j,k,l,uu_h,vv_h,ww_h,rr_h
               uu=uu_h
               vv=vv_h
               ww=ww_h
               rr=rr_h

            else
               print '(a)','readrestart: Attempting to read incompatable turbulence restart file'
               print '(a,4i6)','readrestart: Dimensions in restart file are:',j,k,l
               close(iunit)
               stop
            endif
         close(iunit)
      else
         print '(a)','readrestart: No restart file for inflow turbulence fields available','turbulence'//cit//'.uf'
         stop
      endif
   endif

   if (nturbines > 0) then
      prefix='theta'
      fname =  trim(directory) // trim(prefix) // trim(suffix) // '_' // trim(cit) // trim(ext)
      inquire(file=trim(fname),exist=ex)
      print '(3a)','reading: ',trim(fname)
      if (ex) then
         open(newunit=iunit,file=trim(fname))
            read(iunit,*)theta
         close(iunit)
      else
         print '(a)','readrestart: No restart file for theta avaialble','theta'//cit//'.dat'
         stop
      endif
   endif

   prefix='restart'
   fname =  trim(directory) // trim(prefix) // trim(suffix) // '_' // trim(cit) // trim(ext)
   inquire(file=trim(fname),exist=ex)
   if (ex) then
      print '(3a)','reading: ',trim(fname)
      open(newunit=iunit,file=trim(fname),form="unformatted", status='unknown')
         read(iunit)i,j,k,n,f_h
         f=f_h
      close(iunit)
   else
      print '(3a)','readrestart: restart file does not exist: restart'//cit//'.uf'
      stop
   endif
   deallocate(f_h)
   deallocate(uu_h)
   deallocate(vv_h)
   deallocate(ww_h)
   deallocate(rr_h)
   return
   998 stop 'readrestart: error during read of turbulence restart field'

end subroutine
end module


