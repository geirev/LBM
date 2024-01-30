module m_readrestart
contains
subroutine readrestart(fname,f)
   use mod_dimensions
   character(len=*), intent(in) :: fname
   real, intent(out) :: f(nx,ny,nz,nl)

   logical ex
   integer :: irec,i,j,k,n

   inquire(file=trim(fname),exist=ex)
   if (.not.ex) then
      print '(a)','readrestart: restart file does not exist ',trim(fname)
      stop
   endif

   inquire(iolength=irec)i,j,k,n,f
   open(10,file=trim(fname),form="unformatted", access="direct", recl=irec)
      read(10,rec=1,err=999)i,j,k,n
      if ((i==nx).and.(j==ny).and.(k==nz).and.(n==nl)) then
         read(10,rec=1,err=999)i,j,k,n,f
      else
         print '(a)','readrestart: Attempting to read incompatable restart file'
         print '(a,4i5)','readrestart: Dimensions in restart file are:',i,j,k,nl
         close(10)
         stop
      endif
   close(10)
   return
   999 stop 'readrestart: error during read'

end subroutine
end module


