module m_gnuplotuv
contains
subroutine gnuplotuv(fname,u,v,nx,ny)
   integer, intent(in) :: nx
   integer, intent(in) :: ny
   character(len=*), intent(in) :: fname
   real, intent(in) :: u(nx,ny)
   real, intent(in) :: v(nx,ny)
   integer j,i

   open(10,file=trim(fname))
      do j=1,ny,5
         do i=1,nx,10
            write(10,'(i5,1024g12.4)')i,j,u(i,j),v(i,j)
         enddo
      enddo
   close(10)

end subroutine
end module
