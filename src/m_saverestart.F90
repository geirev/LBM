module m_saverestart
contains
subroutine saverestart(fname,f)
   use mod_dimensions
   character(len=*), intent(in) :: fname
   real, intent(in) :: f(0:nx+1,ny,nz,nl)

   integer :: irec

   inquire(iolength=irec) nx,ny,nz,nl,f
   open(10,file=trim(fname),form="unformatted", access="direct", recl=irec)
      write(10,rec=1)nx,ny,nz,nl,f
   close(10)

end subroutine
end module

