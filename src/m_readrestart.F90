module m_readrestart
contains
subroutine readrestart(f,it,theta)
   use mod_dimensions
   real,    intent(out) :: f(0:nx+1,0:ny+1,0:nz+1,nl)
   integer, intent(in)  :: it
   real,    intent(out) :: theta

   logical ex
   integer :: irec,i,j,k,n
   character(len=6) cit

   write(cit,'(i6.6)')it

   inquire(file='theta'//cit//'.dat',exist=ex)
   if (ex) then
      open(10,file='theta'//cit//'.dat')
         read(10,*)theta
      close(10)
   endif

   inquire(file='restart'//cit//'.uf',exist=ex)
   if (.not.ex) then
      print '(3a)','readrestart: restart file does not exist: restart'//cit//'.uf'
      stop
   endif

   inquire(iolength=irec)i,j,k,n,f
   open(10,file='restart'//cit//'.uf',form="unformatted", access="direct", recl=irec)
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


