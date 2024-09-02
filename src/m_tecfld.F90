module m_tecfld
contains
subroutine tecfld(fname,ja,jb,ka,kb,u,v,w,Ti)
   implicit none
   character(len=*), intent(in) :: fname
   integer,          intent(in) :: ja,jb,ka,kb
   real,             intent(in) :: u(ja:jb,ka:kb)
   real,             intent(in) :: v(ja:jb,ka:kb)
   real,             intent(in) :: w(ja:jb,ka:kb)
   real,             intent(in) :: Ti(ja:jb,ka:kb)
   integer j,k
   logical lopen
   integer iunit
   character(len=100) fn

   do iunit=10,20
      inquire(unit=iunit,opened=lopen)
      if (.not.lopen) then
         open(iunit,file=trim(fname),status='unknown')
            write(iunit,*)'TITLE = "',fname,'"'
            write(iunit,*)'VARIABLES = "i-index" "j-index" "u" "v" "w" "Ti"'
            write(iunit,'(a,i5,a,i5,a)')' ZONE  F=BLOCK, I=',jb-ja+1,', J=',kb-ka+1,', K=1'
            write(iunit,'(30I5)')((j,j=ja,jb),k=ka,kb)
            write(iunit,'(30I5)')((k,j=ja,jb),k=ka,kb)
            write(iunit,900)((u(j,k),j=ja,jb),k=ka,kb)
            write(iunit,900)((v(j,k),j=ja,jb),k=ka,kb)
            write(iunit,900)((w(j,k),j=ja,jb),k=ka,kb)
            write(iunit,900)((Ti(j,k),j=ja,jb),k=ka,kb)
         close(iunit)
         exit
      endif
   enddo
 900 format(10(1x,e12.5))
end subroutine tecfld

end module m_tecfld
