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

   do iunit=10,20
      inquire(unit=iunit,opened=lopen)
      if (.not.lopen) then
         fn(:)=' '
         fn(1:4)='tec_'
         fn(5:5+len_trim(fname))=trim(fname)
         fn(5+len_trim(fname):5+len_trim(fname)+4)='.dat'
         print *,trim(fn)
         open(iunit,file=trim(fn),status='unknown')
            write(iunit,*)'TITLE = "',fname,'"'
            write(iunit,*)'VARIABLES = "i" "j" "k"'
            write(iunit,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj,', K=1'
            write(iunit,'(30I5)')((i,i=1,ii),j=1,jj)
            write(iunit,'(30I5)')((j,i=1,ii),j=1,jj)
            write(iunit,900)((fld(i,j,1),i=1,ii),j=1,jj)
         close(iunit)
         exit
      endif
   enddo
 900 format(10(1x,e12.5))
end subroutine tecfld


end module m_tecfld
