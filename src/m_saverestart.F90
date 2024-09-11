module m_saverestart
contains
subroutine saverestart(f,it,theta)
   use mod_dimensions
   real, intent(in)    :: f(0:nx+1,0:ny+1,0:nz+1,nl)
   integer, intent(in) :: it
   real, intent(in)    :: theta

   character(len=6) cit
   integer :: irec

   write(cit,'(i6.6)')it
   open(10,file='theta'//cit//'.dat')
      write(10,*)theta
   close(10)

   inquire(iolength=irec) nx,ny,nz,nl,f
   open(10,file='restart'//cit//'.uf',form="unformatted", access="direct", recl=irec)
      write(10,rec=1)nx,ny,nz,nl,f
   close(10)

end subroutine
end module

