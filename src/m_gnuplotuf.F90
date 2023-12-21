module m_gnuplotuf
contains
subroutine gnuplotuf(fname,variable,nx,ny)
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   character(len=*), intent(in) :: fname
   real, intent(in) :: variable(nx,ny)
   integer j,i

   integer :: irec

   inquire(iolength=irec) variable
   open(10,file=trim(fname),form="unformatted", access="direct", recl=irec)
      write(10,rec=1)real(variable)
   close(10)

end subroutine
end module
