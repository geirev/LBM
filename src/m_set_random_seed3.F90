module m_set_random_seed3
contains
subroutine set_random_seed3
! Sets a random seed based on the system and wall clock time
   implicit none

   integer , dimension(8)::val
   integer cnt
   integer sze,i
   integer, allocatable, dimension(:):: pt
   logical ex

   call DATE_AND_TIME(values=val)!; print *,'VAL:',val(:)
   call SYSTEM_CLOCK(count=cnt)!; print *,'CNT:',cnt
   call RANDOM_SEED(size=sze)!  ; print *,'sze:',sze
   allocate(pt(sze))
   do i=1,sze,2
      pt(i) = val(8)*val(7)
      pt(i+1) = cnt
   enddo
   inquire(file='seed.dat',exist=ex)
   if (ex) then
      open(10,file='seed.dat')
         read(10,*)pt
      close(10)
   else
      open(10,file='seed.dat')
         write(10,*)pt
      close(10)
   endif
   call RANDOM_SEED(put=pt)
   deallocate(pt)
end subroutine set_random_seed3
end module m_set_random_seed3
