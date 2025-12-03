module m_saverestart
contains
subroutine saverestart(it,f,theta,uu,vv,ww,rr,pottemp,tracer)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile, only : inflowturbulence,nturbines,nrturb,lnodump,iablvisc
#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif
   implicit none
   integer, intent(in)  :: it
   real,    intent(in)  :: theta(nturbines)
   real,    intent(in)  :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in)  :: uu(ny,nz,0:nrturb)
   real,    intent(in)  :: vv(ny,nz,0:nrturb)
   real,    intent(in)  :: ww(ny,nz,0:nrturb)
   real,    intent(in)  :: rr(ny,nz,0:nrturb)
   real,    intent(in)  :: pottemp(:,:,:)
   real,    intent(in)  :: tracer(:,:,:,:)
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

   integer            iunit
   integer ir

   character(len=6)   cit
   character(len=4)   ctile
   character(len=3)   ext
   character(len=10)  prefix
   character(len=10)  directory
   character(len=100) fname

   if (lnodump) return

! File names
   ir=0
#ifdef MPI
   ir=mpi_rank
#endif
   write(ctile,'(i4.4)')ir

   ext='.uf'
   write(cit,'(i6.6)')it
   print '(4a)',' saverestart: tile=',trim(ctile),' iteration=',trim(cit)

   directory='restart/'
   call system('mkdir -p '//trim(directory))

   if (inflowturbulence) then
      prefix='turbulence'
      fname = trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
      open(newunit=iunit,file=trim(fname),form="unformatted", status='replace')
#ifdef _CUDA
         if (.not. allocated(uu_h)) allocate(uu_h(ny,nz,0:nrturb))
         if (.not. allocated(vv_h)) allocate(vv_h(ny,nz,0:nrturb))
         if (.not. allocated(ww_h)) allocate(ww_h(ny,nz,0:nrturb))
         if (.not. allocated(rr_h)) allocate(rr_h(ny,nz,0:nrturb))
         uu_h=uu
         vv_h=vv
         ww_h=ww
         rr_h=rr
         write(iunit)ny,nz,nrturb,uu_h,vv_h,ww_h,rr_h
         deallocate(uu_h,vv_h,ww_h,rr_h)
#else
         write(iunit)ny,nz,nrturb,uu,vv,ww,rr
#endif
      close(iunit)
   endif

   if (ntracer > 0) then
      prefix='tracer'
      fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
      open(newunit=iunit,file=trim(fname),form="unformatted", status='replace')
#ifdef _CUDA
         if (.not. allocated(tracer_h)) allocate(tracer_h(ntracer,0:nx+1,0:ny+1,0:nz+1))
         tracer_h=tracer
         write(iunit)ntracer,tracer_h
         deallocate(tracer_h)
#else
         write(iunit)ntracer,tracer
#endif
      close(iunit)
   endif

   if (iablvisc == 2) then
      prefix='pottemp'
      fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
      open(newunit=iunit,file=trim(fname),form="unformatted", status='replace')
#ifdef _CUDA
         if (.not. allocated(pottemp_h)) allocate(pottemp_h(0:nx+1,0:ny+1,0:nz+1))
         pottemp_h=pottemp
         write(iunit)pottemp_h
         deallocate(pottemp_h)
#else
         write(iunit)pottemp
#endif
      close(iunit)
   endif

   if (nturbines > 0) then
      prefix='theta'
      fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
      open(newunit=iunit,file=trim(fname),form="unformatted", status='replace')
         write(iunit)theta
      close(iunit)
   endif

   prefix='restart'
   fname =  trim(directory) // trim(prefix) // '_' // trim(ctile) // '_' // trim(cit) // trim(ext)
   open(newunit=iunit,file=trim(fname),form="unformatted", status='replace')
#ifdef _CUDA
   if (.not. allocated(f_h)) allocate(f_h(nl,0:nx+1,0:ny+1,0:nz+1))
   f_h = f
   write(iunit) nx,ny,nz,nl,f_h
   deallocate(f_h)
#else
    write(iunit) nx,ny,nz,nl,f
#endif
   close(iunit)

end subroutine
end module

