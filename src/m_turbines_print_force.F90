module m_turbines_print_force
! Dumping detailed diagnostics for one turbine case
contains
subroutine turbines_print_force(force,it)
   use mod_dimensions, only : ny,nz
   use m_turbines_init, only : ieps
   implicit none
   real,    intent(in) :: force(0:ieps,ny,nz,3)
   integer, intent(in) :: it
   character(len=3) tag3
   integer i,j,k

      do i=0,ieps
         write(tag3,'(i3.3)')k
         if (it==1) then
            open(21,file='dist'//tag3//'.dat',status='unknown')
               write(21,*)'TITLE = "Dist',tag3,'"'
               write(21,*)'VARIABLES = "j-index" "k-index" "u" "v" "w" "u2"'
            close(21)
         endif

         open(21,file='dist'//tag3//'.dat',position='append')
            write(21,'(3(a,i5),a)')' ZONE  T="k=',k,'" F=BLOCK, I=',ny,', J=',nz,', K=1'
            write(21,'(30I5)')((j,j=1,ny),k=1,nz)
            write(21,'(30I5)')((k,j=1,ny),k=1,nz)
            write(21,'(10(1x,e12.5))')((force(i,j,k,1),j=1,ny),k=1,nz)
            write(21,'(10(1x,e12.5))')((force(i,j,k,2),j=1,ny),k=1,nz)
            write(21,'(10(1x,e12.5))')((force(i,j,k,3),j=1,ny),k=1,nz)
            write(21,'(10(1x,e12.5))')((sqrt(force(i,j,k,2)**2+force(i,j,k,3)**2),j=1,ny),k=1,nz)
         close(21)
      enddo
end subroutine
end module

