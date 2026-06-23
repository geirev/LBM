module m_tecfld
contains
subroutine tecfld(fname,vname,ii,jj,nr,fld,ir,np)
   implicit none
   character(len=*), intent(in) :: fname
   character(len=*), intent(in) :: vname
   integer,          intent(in) :: ii,jj,nr,ir,np
   real,             intent(in) :: fld(ii,jj+1,nr)
   integer i,j,ja,jb,jj1
   logical lopen
   integer iunit
   character(len=100) fn
   fn=trim(fname)//'.dat'
   print *,trim(fn)
   open(newunit=iunit,file=trim(fn),status='unknown')
      ja=ir*jj+1
      jb=ir*jj+jj
      jj1=jj
      if (ir < np - 1) then
         jb=jb+1
         jj1=jj+1
      endif
      write(iunit,*)'TITLE = "',fname,'"'
      write(iunit,*)'VARIABLES = "i" "j" "',trim(vname),'"'
      write(iunit,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',ii,', J=',jj1,', K=1'
      write(iunit,'(30I5)')((i,i=1,ii),j=ja,jb)
      write(iunit,'(30I5)')((j,i=1,ii),j=ja,jb)
      write(iunit,900)((fld(i,j,1),i=1,ii),j=1,jj1)
   close(iunit)
 900 format(10(1x,e12.5))
end subroutine tecfld
end module m_tecfld
