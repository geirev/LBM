module m_tecfld
contains
subroutine tecfld(fname,ii,jj,nr,fld)
   implicit none
   character(len=*), intent(in) :: fname
   integer,          intent(in) :: ii,jj,nr
   real,             intent(in) :: fld(ii,jj,nr)
   integer i,j
   logical lopen
   integer iunit
   character(len=100) fn
   fn='tec_'//trim(fname)//'.dat'
   print *,trim(fn)
   open(newunit=iunit,file=trim(fn),status='unknown')
      write(iunit,*)'TITLE = "',fname,'"'
      write(iunit,*)'VARIABLES = "i" "j" "',trim(fname),'"'
      write(iunit,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'
      write(iunit,'(30I5)')((i,i=1,ii),j=1,jj)
      write(iunit,'(30I5)')((j,i=1,ii),j=1,jj)
      write(iunit,900)((fld(i,j,1),i=1,ii),j=1,jj)
   close(iunit)
 900 format(10(1x,e12.5))
end subroutine tecfld


end module m_tecfld
