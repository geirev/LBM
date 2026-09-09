module m_saverestart
   implicit none
   private
   public :: saverestart

contains

subroutine ioserror(ios,ir,message,fname,iomsg,iunit)

   implicit none

   integer,          intent(in) :: ios
   integer,          intent(in) :: ir
   character(len=*), intent(in) :: message
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: iomsg
   integer, optional,intent(in) :: iunit

   if (ios == 0) return

   print *
   print *, '============================================================'
   print *, 'saverestart: I/O error'
   print *, 'operation : ', trim(message)
   print *, 'rank      : ', ir
   print *, 'file      : ', trim(fname)
   print *, 'iostat    : ', ios
   print *, 'message   : ', trim(iomsg)
   print *, '============================================================'
   print *

   if (present(iunit)) close(iunit)

   error stop

end subroutine ioserror


subroutine allocerror(istat,ir,message,errmsg)

   implicit none

   integer,          intent(in) :: istat
   integer,          intent(in) :: ir
   character(len=*), intent(in) :: message
   character(len=*), intent(in) :: errmsg

   if (istat == 0) return

   print *
   print *, '============================================================'
   print *, 'saverestart: allocation error'
   print *, 'operation : ', trim(message)
   print *, 'rank      : ', ir
   print *, 'stat      : ', istat
   print *, 'message   : ', trim(errmsg)
   print *, '============================================================'
   print *

   error stop

end subroutine allocerror


subroutine saverestart(it,f,uu,vv,ww,rr,pottemp,tracer)

   use mod_dimensions
   use mod_turbines,   only : turbines
   use mod_D3Q27setup, only : nl
   use m_readinfile,   only : inflowturbulence, nturbines, nrturb, &
                              ldump, iablvisc

#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif

   implicit none

   integer, intent(in) :: it

   real, intent(in) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in) :: uu(ny,nz,0:nrturb)
   real, intent(in) :: vv(ny,nz,0:nrturb)
   real, intent(in) :: ww(ny,nz,0:nrturb)
   real, intent(in) :: rr(ny,nz,0:nrturb)
   real, intent(in) :: pottemp(:,:,:)
   real, intent(in) :: tracer(:,:,:,:)

#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: uu
   attributes(device) :: vv
   attributes(device) :: ww
   attributes(device) :: rr
   attributes(device) :: tracer
   attributes(device) :: pottemp

   real, allocatable :: f_h(:,:,:,:)
   real, allocatable :: uu_h(:,:,:)
   real, allocatable :: vv_h(:,:,:)
   real, allocatable :: ww_h(:,:,:)
   real, allocatable :: rr_h(:,:,:)
   real, allocatable :: tracer_h(:,:,:,:)
   real, allocatable :: pottemp_h(:,:,:)
#endif

   integer :: iunit
   integer :: ir
   integer :: ios
   integer :: istat

   character(len=6)   :: cit
   character(len=4)   :: ctile
   character(len=3)   :: ext
   character(len=10)  :: prefix
   character(len=10)  :: directory
   character(len=100) :: fname
   character(len=256) :: iomsg
   character(len=256) :: errmsg


!-----------------------------------------------------------------------
! Restart output enabled?
!-----------------------------------------------------------------------

   if (.not. ldump) return


!-----------------------------------------------------------------------
! MPI tile and iteration
!-----------------------------------------------------------------------

   ir = 0

#ifdef MPI
   ir = mpi_rank
#endif

   write(ctile,'(i4.4)') ir

   ext = '.uf'
   write(cit,'(i6.6)') it

   print '(4a)', ' saverestart: tile=', trim(ctile), &
                  ' iteration=', trim(cit)

   directory = 'restart/'

   call system('mkdir -p '//trim(directory))


!=======================================================================
! Inflow turbulence
!
! record 1 : ny, nz, nrturb
! record 2 : uu, vv, ww, rr
!=======================================================================

   if (inflowturbulence) then

      prefix = 'turbulence'

      fname = trim(directory) // trim(prefix) // '_' // &
              trim(ctile) // '_' // trim(cit) // trim(ext)

      print '(3a)', 'writing: ', trim(fname)

      open(newunit=iunit, file=trim(fname), form='unformatted', &
           status='replace', action='write', &
           iostat=ios, iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'opening turbulence restart',fname,iomsg)

      write(iunit,iostat=ios,iomsg=iomsg) ny,nz,nrturb

      if (ios /= 0) &
         call ioserror(ios,ir,'writing turbulence header', &
                       fname,iomsg,iunit)

#ifdef _CUDA

      allocate(uu_h(ny,nz,0:nrturb),stat=istat,errmsg=errmsg)
      if (istat /= 0) &
         call allocerror(istat,ir,'allocating uu_h',errmsg)

      allocate(vv_h(ny,nz,0:nrturb),stat=istat,errmsg=errmsg)
      if (istat /= 0) &
         call allocerror(istat,ir,'allocating vv_h',errmsg)

      allocate(ww_h(ny,nz,0:nrturb),stat=istat,errmsg=errmsg)
      if (istat /= 0) &
         call allocerror(istat,ir,'allocating ww_h',errmsg)

      allocate(rr_h(ny,nz,0:nrturb),stat=istat,errmsg=errmsg)
      if (istat /= 0) &
         call allocerror(istat,ir,'allocating rr_h',errmsg)

      uu_h = uu
      vv_h = vv
      ww_h = ww
      rr_h = rr

      write(iunit,iostat=ios,iomsg=iomsg) uu_h,vv_h,ww_h,rr_h

      if (ios /= 0) &
         call ioserror(ios,ir,'writing turbulence fields', &
                       fname,iomsg,iunit)

      deallocate(uu_h,vv_h,ww_h,rr_h)

#else

      write(iunit,iostat=ios,iomsg=iomsg) uu,vv,ww,rr

      if (ios /= 0) &
         call ioserror(ios,ir,'writing turbulence fields', &
                       fname,iomsg,iunit)

#endif

      close(iunit,iostat=ios,iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'closing turbulence restart',fname,iomsg)

   endif


!=======================================================================
! Tracer
!
! record 1 : ntracer
! record 2 : tracer
!=======================================================================

   if (ntracer > 0) then

      prefix = 'tracer'

      fname = trim(directory) // trim(prefix) // '_' // &
              trim(ctile) // '_' // trim(cit) // trim(ext)

      print '(3a)', 'writing: ', trim(fname)

      open(newunit=iunit, file=trim(fname), form='unformatted', &
           status='replace', action='write', &
           iostat=ios, iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'opening tracer restart',fname,iomsg)

      write(iunit,iostat=ios,iomsg=iomsg) ntracer

      if (ios /= 0) &
         call ioserror(ios,ir,'writing tracer header', &
                       fname,iomsg,iunit)

#ifdef _CUDA

      allocate(tracer_h(ntracer,0:nx+1,0:ny+1,0:nz+1), &
               stat=istat,errmsg=errmsg)

      if (istat /= 0) &
         call allocerror(istat,ir,'allocating tracer_h',errmsg)

      tracer_h = tracer

      write(iunit,iostat=ios,iomsg=iomsg) tracer_h

      if (ios /= 0) &
         call ioserror(ios,ir,'writing tracer fields', &
                       fname,iomsg,iunit)

      deallocate(tracer_h)

#else

      write(iunit,iostat=ios,iomsg=iomsg) tracer

      if (ios /= 0) &
         call ioserror(ios,ir,'writing tracer fields', &
                       fname,iomsg,iunit)

#endif

      close(iunit,iostat=ios,iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'closing tracer restart',fname,iomsg)

   endif


!=======================================================================
! Potential temperature
!
! record 1 : nx, ny, nz
! record 2 : pottemp
!=======================================================================

   if (iablvisc == 2) then

      prefix = 'pottemp'

      fname = trim(directory) // trim(prefix) // '_' // &
              trim(ctile) // '_' // trim(cit) // trim(ext)

      print '(3a)', 'writing: ', trim(fname)

      open(newunit=iunit, file=trim(fname), form='unformatted', &
           status='replace', action='write', &
           iostat=ios, iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'opening pottemp restart',fname,iomsg)

      write(iunit,iostat=ios,iomsg=iomsg) nx,ny,nz

      if (ios /= 0) &
         call ioserror(ios,ir,'writing pottemp header', &
                       fname,iomsg,iunit)

#ifdef _CUDA

      allocate(pottemp_h(0:nx+1,0:ny+1,0:nz+1), &
               stat=istat,errmsg=errmsg)

      if (istat /= 0) &
         call allocerror(istat,ir,'allocating pottemp_h',errmsg)

      pottemp_h = pottemp

      write(iunit,iostat=ios,iomsg=iomsg) pottemp_h

      if (ios /= 0) &
         call ioserror(ios,ir,'writing potential temperature', &
                       fname,iomsg,iunit)

      deallocate(pottemp_h)

#else

      write(iunit,iostat=ios,iomsg=iomsg) pottemp

      if (ios /= 0) &
         call ioserror(ios,ir,'writing potential temperature', &
                       fname,iomsg,iunit)

#endif

      close(iunit,iostat=ios,iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'closing pottemp restart',fname,iomsg)

   endif


!=======================================================================
! Turbines
!
! record 1 : nturbines
! record 2 : turbines
!=======================================================================

   if (nturbines > 0) then

      prefix = 'turbines'

      fname = trim(directory) // trim(prefix) // '_' // &
              trim(ctile) // '_' // trim(cit) // trim(ext)

      print '(3a)', 'writing: ', trim(fname)

      open(newunit=iunit, file=trim(fname), form='unformatted', &
           status='replace', action='write', &
           iostat=ios, iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'opening turbine restart',fname,iomsg)

      write(iunit,iostat=ios,iomsg=iomsg) nturbines

      if (ios /= 0) &
         call ioserror(ios,ir,'writing turbine header', &
                       fname,iomsg,iunit)

      write(iunit,iostat=ios,iomsg=iomsg) turbines

      if (ios /= 0) &
         call ioserror(ios,ir,'writing turbine data', &
                       fname,iomsg,iunit)

      close(iunit,iostat=ios,iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'closing turbine restart',fname,iomsg)

   endif


!=======================================================================
! Main LBM restart
!
! record 1 : nx, ny, nz, nl
! record 2 : f
!=======================================================================

   prefix = 'restart'

   fname = trim(directory) // trim(prefix) // '_' // &
           trim(ctile) // '_' // trim(cit) // trim(ext)

   print '(3a)', 'writing: ', trim(fname)

   open(newunit=iunit, file=trim(fname), form='unformatted', &
        status='replace', action='write', &
        iostat=ios, iomsg=iomsg)

   if (ios /= 0) &
      call ioserror(ios,ir,'opening restart file',fname,iomsg)

   write(iunit,iostat=ios,iomsg=iomsg) nx,ny,nz,nl

   if (ios /= 0) &
      call ioserror(ios,ir,'writing restart header', &
                    fname,iomsg,iunit)

#ifdef _CUDA

   allocate(f_h(nl,0:nx+1,0:ny+1,0:nz+1), &
            stat=istat,errmsg=errmsg)

   if (istat /= 0) &
      call allocerror(istat,ir,'allocating f_h',errmsg)

   f_h = f

   write(iunit,iostat=ios,iomsg=iomsg) f_h

   if (ios /= 0) &
      call ioserror(ios,ir,'writing distribution functions', &
                    fname,iomsg,iunit)

   deallocate(f_h)

#else

   write(iunit,iostat=ios,iomsg=iomsg) f

   if (ios /= 0) &
      call ioserror(ios,ir,'writing distribution functions', &
                    fname,iomsg,iunit)

#endif

   close(iunit,iostat=ios,iomsg=iomsg)

   if (ios /= 0) &
      call ioserror(ios,ir,'closing restart file',fname,iomsg)

   print *, 'saverestart: completed, rank=',ir


end subroutine saverestart

end module m_saverestart
