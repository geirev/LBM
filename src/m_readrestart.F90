module m_readrestart
contains
subroutine readrestart(it,f,theta,uu,vv,ww,rr,pottemp,tracer)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile, only : inflowturbulence,nturbines,nrturb,iablvisc
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif
   implicit none
   integer, intent(in)  :: it
   real,    intent(out) :: theta(nturbines)
   real,    intent(out) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(out) :: uu(ny,nz,0:nrturb)
   real,    intent(out) :: vv(ny,nz,0:nrturb)
   real,    intent(out) :: ww(ny,nz,0:nrturb)
   real,    intent(out) :: rr(ny,nz,0:nrturb)
   real,    intent(out) :: tracer(:,:,:,:)
   real,    intent(out) :: pottemp(:,:,:)
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
   attributes(device) :: rr
   attributes(device) :: tracer
   attributes(device) :: pottemp
#endif

#ifdef _CUDA
   real, allocatable :: f_h(:,:,:,:)
   real, allocatable :: uu_h(:,:,:)
   real, allocatable :: vv_h(:,:,:)
   real, allocatable :: ww_h(:,:,:)
   real, allocatable :: rr_h(:,:,:)
   real, allocatable :: tracer_h(:,:,:,:)
   real, allocatable :: pottemp_h(:,:,:)
#endif

   logical ex
   integer :: i,j,k,l,n
   integer iunit
   integer ir

   character(len=6) cit
   character(len=4) ctile
   character(len=3) ext
   character(len=10) prefix
   character(len=10)  directory
   character(len=100) fname

! File names
   ir=0
#ifdef MPI
   ir=mpi_rank
#endif
   write(ctile,'(i4.4)')ir

   ext='.uf'
   write(cit,'(i6.6)')it
   print '(4a)',' readrestart: tile=',trim(ctile),' iteration=',trim(cit)

   directory='restart/'
   call system('mkdir -p '//trim(directory))

   if (inflowturbulence) then
      prefix='turbulence'
      fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
      inquire(file=trim(fname),exist=ex)
      print '(3a)','reading: ',trim(fname)
      if (ex) then
         open(newunit=iunit,file=trim(fname),form="unformatted", status='old')
            read(iunit,err=998)j,k,l
            rewind(iunit)
            if ((j==ny).and.(k==nz).and.(l==nrturb)) then
#ifdef _CUDA
               if (.not. allocated(uu_h)) allocate(uu_h(ny,nz,0:nrturb))
               if (.not. allocated(vv_h)) allocate(vv_h(ny,nz,0:nrturb))
               if (.not. allocated(ww_h)) allocate(ww_h(ny,nz,0:nrturb))
               if (.not. allocated(rr_h)) allocate(rr_h(ny,nz,0:nrturb))
               read(iunit,err=998)j,k,l,uu_h,vv_h,ww_h,rr_h
               uu=uu_h
               vv=vv_h
               ww=ww_h
               rr=rr_h
               deallocate(uu_h,vv_h,ww_h,rr_h)
#else
               read(iunit,err=998)j,k,l,uu,vv,ww,rr
#endif
            else
               print '(a)','readrestart: Attempting to read incompatable turbulence restart file'
               print '(a,4i6)','readrestart: Dimensions in restart file are:',j,k,l
               close(iunit)
               stop
            endif
         close(iunit)
      else
         print '(a)','readrestart: No restart file for inflow turbulence fields available',trim(fname)
         stop
      endif
   endif

   if (nturbines > 0) then
      prefix='theta'
      fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
      inquire(file=trim(fname),exist=ex)
      print '(3a)','reading: ',trim(fname)
      if (ex) then
         open(newunit=iunit,file=trim(fname),form="unformatted", status='unknown')
            read(iunit,err=998)theta
         close(iunit)
         print *,'read restart theta',theta
      else
         print '(a)','readrestart: No restart file for theta avaialble:',trim(fname)
         stop
      endif
   endif

   if (iablvisc == 2) then
      prefix='pottemp'
      fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
      inquire(file=trim(fname),exist=ex)
      print '(3a)','reading: ',trim(fname)
      if (ex) then
         open(newunit=iunit,file=trim(fname),form="unformatted", status='unknown')
#ifdef _CUDA
         if (.not. allocated(pottemp_h)) allocate(pottemp_h(0:nx+1,0:ny+1,0:nz+1))
         read(iunit)pottemp_h
         pottemp=pottemp_h
         deallocate(pottemp_h)
#else
         read(iunit)pottemp
#endif
         close(iunit)
      else
         print '(a)','readrestart: No restart file for potential temperature avaialble:',trim(fname)
         stop
      endif
   endif

   if (ntracer > 0) then
      prefix='tracer'
      fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
      inquire(file=trim(fname),exist=ex)
      print '(3a)','reading: ',trim(fname)
      if (ex) then
         open(newunit=iunit,file=trim(fname),form="unformatted", status='unknown')
#ifdef _CUDA
         if (.not. allocated(tracer_h)) allocate(tracer_h(ntracer,0:nx+1,0:ny+1,0:nz+1))
         read(iunit)n,tracer_h
         tracer=tracer_h
         deallocate(tracer_h)
#else
         read(iunit)n,tracer
#endif
         close(iunit)
      else
         print '(a)','readrestart: No restart file for tracer avaialble:',trim(fname)
         stop
      endif
   endif

   prefix='restart'
   fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
   inquire(file=trim(fname),exist=ex)
   if (ex) then
      print '(3a)','reading: ',trim(fname)
      open(newunit=iunit,file=trim(fname),form="unformatted", status='unknown')
         read(iunit)i,j,k,l
         rewind(iunit)
         if ((i==nx).and.(j==ny).and.(k==nz).and.(l==nl)) then
#ifdef _CUDA
            allocate(f_h(nl,0:nx+1,0:ny+1,0:nz+1))
            read(iunit,err=998)i,j,k,l,f_h
            f=f_h
            deallocate(f_h)
#else
            read(iunit,err=998)i,j,k,l,f
#endif
         endif
      close(iunit)
   else
      print '(3a)','readrestart: restart file does not exist: ',trim(fname)
      stop
   endif
   return
   998 stop 'readrestart: error during read of turbulence restart field'

end subroutine
end module


