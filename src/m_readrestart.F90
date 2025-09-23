module m_readrestart
contains
subroutine readrestart(it,f,theta,uu,vv,ww,rr)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_readinfile, only : inflowturbulence,nturbines,nrturb
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

   allocate(f_h(nl,0:nx+1,0:ny+1,0:nz+1))
   allocate(uu_h(ny,nz,0:nrturb)        )
   allocate(vv_h(ny,nz,0:nrturb)        )
   allocate(ww_h(ny,nz,0:nrturb)        )
   allocate(rr_h(ny,nz,0:nrturb)        )

   write(cit,'(i6.6)')it
   print '(a,a)','readrestart:',cit

   if (inflowturbulence) then
      print '(3a)','reading: turbulence'//cit//'.uf'
      inquire(file='turbulence'//cit//'.uf',exist=ex)
      if (ex) then
         open(newunit=iunit,file='turbulence'//cit//'.uf',form="unformatted", status='old')
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
      print '(3a)','reading: theta'//cit//'.dat'
      inquire(file='theta'//cit//'.dat',exist=ex)
      if (ex) then
         open(newunit=iunit,file='theta'//cit//'.dat')
            read(iunit,*)theta
         close(iunit)
      else
         print '(a)','readrestart: No restart file for theta avaialble','theta'//cit//'.dat'
         stop
      endif
   endif

   inquire(file='restart'//cit//'.uf',exist=ex)
   if (ex) then
      open(newunit=iunit,file='restart'//cit//'.uf',form="unformatted", status='unknown')
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


