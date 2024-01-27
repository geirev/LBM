module m_gnuplot
contains
subroutine gnuplot(fname,variable,nx,ny)
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   character(len=*), intent(in) :: fname
   real, intent(in) :: variable(nx,ny)
   integer j,i

   integer :: irec

   open(10,file=trim(fname))
      write(10,'(1025i5)')nx,(i,i=1,nx)
      do j=1,ny
         write(10,'(i5,1024g12.4)')j,variable(:,j)
      enddo
   close(10)

end subroutine
end module
