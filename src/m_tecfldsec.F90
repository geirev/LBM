module m_tecfldsec
contains
subroutine tecfldsec(fname,nrsec,ja,jb,ka,kb,u,v,w,Ti)
   implicit none
   character(len=*), intent(in) :: fname
   integer,          intent(in) :: nrsec
   integer,          intent(in) :: ja,jb,ka,kb
   real,             intent(in) :: u(0:nrsec-1,ja:jb,ka:kb)
   real,             intent(in) :: v(0:nrsec-1,ja:jb,ka:kb)
   real,             intent(in) :: w(0:nrsec-1,ja:jb,ka:kb)
   real,             intent(in) :: Ti(0:nrsec-1,ja:jb,ka:kb)
   integer j,k,isec
   logical lopen
   integer iunit
   character(len=2) csec

   do iunit=10,20
      inquire(unit=iunit,opened=lopen)
      if (.not.lopen) then
         open(iunit,file=trim(fname),status='unknown')
            write(iunit,*)'TITLE = "',fname,'"'
            write(iunit,*)'VARIABLES = "j-index" "k-index" "u" "v" "w" "Ti"'
            do isec=0,nrsec-1
               write(csec,'(i2.2)')isec
               if (isec==0) then
                  write(iunit,'(3a,i5,a,i5,a)')' ZONE T="D',csec,'" F=BLOCK, I=',jb-ja+1,', J=',kb-ka+1,', K=1'
                  write(iunit,'(30I5)')((j,j=ja,jb),k=ka,kb)
                  write(iunit,'(30I5)')((k,j=ja,jb),k=ka,kb)
               else
                  write(iunit,'(3a,i5,a,i5,a)')' ZONE T="D',csec,'" F=BLOCK, I=',jb-ja+1,', J=',kb-ka+1,&
                                               &', K=1 VARSHARELIST=([1,2]=1)'
               endif
               write(iunit,900)(( u(isec,j,k),j=ja,jb),k=ka,kb)
               write(iunit,900)(( v(isec,j,k),j=ja,jb),k=ka,kb)
               write(iunit,900)(( w(isec,j,k),j=ja,jb),k=ka,kb)
               write(iunit,900)((Ti(isec,j,k),j=ja,jb),k=ka,kb)
            enddo
         close(iunit)
         exit
      endif
   enddo
 900 format(10(1x,e12.5))
end subroutine

end module
