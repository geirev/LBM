module m_readrestart
contains
subroutine readrestart(it,f,theta,uu,vv,ww,rr)
   use mod_dimensions
   use m_readinfile, only : inflowturbulence,nturbines
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
   integer :: irec,i,j,k,l,n
   character(len=6) cit

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
         inquire(iolength=irec)j,k,l,uu_h,vv_h,ww_h,rr_h
         open(10,file='turbulence'//cit//'.uf',form="unformatted", access="direct", recl=irec)
            read(10,rec=1,err=998)j,k,l
            if ((j==ny).and.(k==nz).and.(l==nrturb)) then
               read(10,rec=1,err=998)j,k,l,uu_h,vv_h,ww_h,rr_h
               uu=uu_h
               vv=vv_h
               ww=ww_h
               rr=rr_h

            else
               print '(a)','readrestart: Attempting to read incompatable turbulence restart file'
               print '(a,4i6)','readrestart: Dimensions in restart file are:',j,k,l
               close(10)
               stop
            endif
         close(10)
      else
         print '(a)','readrestart: No restart file for inflow turbulence fields available','turbulence'//cit//'.uf'
         stop
      endif
   endif

   if (nturbines > 0) then
      print '(3a)','reading: theta'//cit//'.uf'
      inquire(file='theta'//cit//'.dat',exist=ex)
      if (ex) then
         open(10,file='theta'//cit//'.dat')
            read(10,*)theta
         close(10)
      else
         print '(a)','readrestart: No restart file for theta avaialble','theta'//cit//'.dat'
         stop
      endif
   endif

   inquire(file='restart'//cit//'.uf',exist=ex)
   if (ex) then
      inquire(iolength=irec)i,j,k,n,f_h
      open(10,file='restart'//cit//'.uf',form="unformatted", access="direct", recl=irec)
         read(10,rec=1,err=999)i,j,k,n
         if ((i==nx).and.(j==ny).and.(k==nz).and.(n==nl)) then
            read(10,rec=1,err=999)i,j,k,n,f_h
            f=f_h
         else
            print '(a)','readrestart: Attempting to read incompatable restart file'
            print '(a,4i5)','readrestart: Dimensions in restart file are:',i,j,k,nl
            close(10)
            stop
         endif
      close(10)
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
   999 stop 'readrestart: error during read of restart file'

end subroutine
end module


