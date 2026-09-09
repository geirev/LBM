module m_readrestart
   implicit none
   private
   public :: readrestart

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
   print *, 'readrestart: I/O error'
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
   print *, 'readrestart: allocation error'
   print *, 'operation : ', trim(message)
   print *, 'rank      : ', ir
   print *, 'stat      : ', istat
   print *, 'message   : ', trim(errmsg)
   print *, '============================================================'
   print *

   error stop

end subroutine allocerror


subroutine readrestart(it,f,uu,vv,ww,rr,pottemp,tracer)

   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile, only : inflowturbulence, nturbines, nrturb, iablvisc
   use mod_turbines, only : turbines

#ifdef MPI
   use m_mpi_decomp_init, only : mpi_rank
#endif

   implicit none

   integer, intent(in) :: it

   real, intent(out) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(out) :: uu(ny,nz,0:nrturb)
   real, intent(out) :: vv(ny,nz,0:nrturb)
   real, intent(out) :: ww(ny,nz,0:nrturb)
   real, intent(out) :: rr(ny,nz,0:nrturb)
   real, intent(out) :: tracer(:,:,:,:)
   real, intent(out) :: pottemp(:,:,:)

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

   logical :: ex

   integer :: i,j,k,l,n
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
! MPI tile and iteration
!-----------------------------------------------------------------------

   ir = 0

#ifdef MPI
   ir = mpi_rank
#endif

   write(ctile,'(i4.4)') ir

   ext = '.uf'
   write(cit,'(i6.6)') it

   print '(4a)', ' readrestart: tile=', trim(ctile), &
                  ' iteration=', trim(cit)

   directory = 'restart/'


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

      print '(3a)', 'reading: ', trim(fname)

      inquire(file=trim(fname),exist=ex)

      if (.not. ex) then
         print *, 'readrestart: turbulence restart file does not exist'
         print *, 'rank = ', ir
         print *, 'file = ', trim(fname)
         error stop
      endif

      open(newunit=iunit,file=trim(fname),form='unformatted', &
           status='old',action='read',iostat=ios,iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'opening turbulence restart',fname,iomsg)

      read(iunit,iostat=ios,iomsg=iomsg) j,k,l

      if (ios /= 0) &
         call ioserror(ios,ir,'reading turbulence header', &
                       fname,iomsg,iunit)

      if ((j /= ny) .or. (k /= nz) .or. (l /= nrturb)) then
         print *
         print *, '============================================================'
         print *, 'readrestart: incompatible turbulence restart dimensions'
         print *, 'rank      : ', ir
         print *, 'file      : ', trim(fname)
         print *, 'file      : ', j,k,l
         print *, 'current   : ', ny,nz,nrturb
         print *, '============================================================'
         print *
         close(iunit)
         error stop
      endif

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

      read(iunit,iostat=ios,iomsg=iomsg) uu_h,vv_h,ww_h,rr_h

      if (ios /= 0) &
         call ioserror(ios,ir,'reading turbulence fields', &
                       fname,iomsg,iunit)

      uu = uu_h
      vv = vv_h
      ww = ww_h
      rr = rr_h

      deallocate(uu_h,vv_h,ww_h,rr_h)

#else

      read(iunit,iostat=ios,iomsg=iomsg) uu,vv,ww,rr

      if (ios /= 0) &
         call ioserror(ios,ir,'reading turbulence fields', &
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

      print '(3a)', 'reading: ', trim(fname)

      inquire(file=trim(fname),exist=ex)

      if (.not. ex) then
         print *, 'readrestart: tracer restart file does not exist'
         print *, 'rank = ', ir
         print *, 'file = ', trim(fname)
         error stop
      endif

      open(newunit=iunit,file=trim(fname),form='unformatted', &
           status='old',action='read',iostat=ios,iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'opening tracer restart',fname,iomsg)

      read(iunit,iostat=ios,iomsg=iomsg) n

      if (ios /= 0) &
         call ioserror(ios,ir,'reading tracer header', &
                       fname,iomsg,iunit)

      if (n /= ntracer) then
         print *
         print *, '============================================================'
         print *, 'readrestart: incompatible tracer restart'
         print *, 'rank             : ', ir
         print *, 'file             : ', trim(fname)
         print *, 'file ntracer     : ', n
         print *, 'current ntracer  : ', ntracer
         print *, '============================================================'
         print *
         close(iunit)
         error stop
      endif

#ifdef _CUDA

      allocate(tracer_h(ntracer,0:nx+1,0:ny+1,0:nz+1), &
               stat=istat,errmsg=errmsg)

      if (istat /= 0) &
         call allocerror(istat,ir,'allocating tracer_h',errmsg)

      read(iunit,iostat=ios,iomsg=iomsg) tracer_h

      if (ios /= 0) &
         call ioserror(ios,ir,'reading tracer fields', &
                       fname,iomsg,iunit)

      tracer = tracer_h

      deallocate(tracer_h)

#else

      read(iunit,iostat=ios,iomsg=iomsg) tracer

      if (ios /= 0) &
         call ioserror(ios,ir,'reading tracer fields', &
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

      print '(3a)', 'reading: ', trim(fname)

      inquire(file=trim(fname),exist=ex)

      if (.not. ex) then
         print *, 'readrestart: pottemp restart file does not exist'
         print *, 'rank = ', ir
         print *, 'file = ', trim(fname)
         error stop
      endif

      open(newunit=iunit,file=trim(fname),form='unformatted', &
           status='old',action='read',iostat=ios,iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'opening pottemp restart',fname,iomsg)

      read(iunit,iostat=ios,iomsg=iomsg) i,j,k

      if (ios /= 0) &
         call ioserror(ios,ir,'reading pottemp header', &
                       fname,iomsg,iunit)

      if ((i /= nx) .or. (j /= ny) .or. (k /= nz)) then
         print *
         print *, '============================================================'
         print *, 'readrestart: incompatible pottemp restart dimensions'
         print *, 'rank      : ', ir
         print *, 'file      : ', trim(fname)
         print *, 'file      : ', i,j,k
         print *, 'current   : ', nx,ny,nz
         print *, '============================================================'
         print *
         close(iunit)
         error stop
      endif

#ifdef _CUDA

      allocate(pottemp_h(0:nx+1,0:ny+1,0:nz+1), &
               stat=istat,errmsg=errmsg)

      if (istat /= 0) &
         call allocerror(istat,ir,'allocating pottemp_h',errmsg)

      read(iunit,iostat=ios,iomsg=iomsg) pottemp_h

      if (ios /= 0) &
         call ioserror(ios,ir,'reading potential temperature', &
                       fname,iomsg,iunit)

      pottemp = pottemp_h

      deallocate(pottemp_h)

#else

      read(iunit,iostat=ios,iomsg=iomsg) pottemp

      if (ios /= 0) &
         call ioserror(ios,ir,'reading potential temperature', &
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

      print '(3a)', 'reading: ', trim(fname)

      inquire(file=trim(fname),exist=ex)

      if (.not. ex) then
         print *, 'readrestart: turbine restart file does not exist'
         print *, 'rank = ', ir
         print *, 'file = ', trim(fname)
         error stop
      endif

      open(newunit=iunit,file=trim(fname),form='unformatted', &
           status='old',action='read',iostat=ios,iomsg=iomsg)

      if (ios /= 0) &
         call ioserror(ios,ir,'opening turbine restart',fname,iomsg)

      read(iunit,iostat=ios,iomsg=iomsg) n

      if (ios /= 0) &
         call ioserror(ios,ir,'reading turbine header', &
                       fname,iomsg,iunit)

      if (n /= nturbines) then
         print *
         print *, '============================================================'
         print *, 'readrestart: incompatible turbine restart'
         print *, 'rank               : ', ir
         print *, 'file               : ', trim(fname)
         print *, 'file nturbines     : ', n
         print *, 'current nturbines  : ', nturbines
         print *, '============================================================'
         print *
         close(iunit)
         error stop
      endif

      read(iunit,iostat=ios,iomsg=iomsg) turbines

      if (ios /= 0) &
         call ioserror(ios,ir,'reading turbine data', &
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

   print '(3a)', 'reading: ', trim(fname)

   inquire(file=trim(fname),exist=ex)

   if (.not. ex) then
      print *, 'readrestart: restart file does not exist'
      print *, 'rank = ', ir
      print *, 'file = ', trim(fname)
      error stop
   endif

   open(newunit=iunit,file=trim(fname),form='unformatted', &
        status='old',action='read',iostat=ios,iomsg=iomsg)

   if (ios /= 0) &
      call ioserror(ios,ir,'opening restart file',fname,iomsg)

   read(iunit,iostat=ios,iomsg=iomsg) i,j,k,l

   if (ios /= 0) &
      call ioserror(ios,ir,'reading restart header', &
                    fname,iomsg,iunit)

   if ((i /= nx) .or. (j /= ny) .or. &
       (k /= nz) .or. (l /= nl)) then

      print *
      print *, '============================================================'
      print *, 'readrestart: incompatible restart dimensions'
      print *, 'rank      : ', ir
      print *, 'file      : ', trim(fname)
      print *, 'file      : ', i,j,k,l
      print *, 'current   : ', nx,ny,nz,nl
      print *, '============================================================'
      print *

      close(iunit)
      error stop

   endif


#ifdef _CUDA

   allocate(f_h(nl,0:nx+1,0:ny+1,0:nz+1), &
            stat=istat,errmsg=errmsg)

   if (istat /= 0) &
      call allocerror(istat,ir,'allocating f_h',errmsg)

   read(iunit,iostat=ios,iomsg=iomsg) f_h

   if (ios /= 0) &
      call ioserror(ios,ir,'reading distribution functions', &
                    fname,iomsg,iunit)

   f = f_h

   deallocate(f_h)

#else

   read(iunit,iostat=ios,iomsg=iomsg) f

   if (ios /= 0) &
      call ioserror(ios,ir,'reading distribution functions', &
                    fname,iomsg,iunit)

#endif


   close(iunit,iostat=ios,iomsg=iomsg)

   if (ios /= 0) &
      call ioserror(ios,ir,'closing restart file',fname,iomsg)

   print *, 'readrestart: completed, rank=',ir


end subroutine readrestart

end module m_readrestart
