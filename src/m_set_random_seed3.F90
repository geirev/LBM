module m_set_random_seed3
contains
subroutine set_random_seed3
! Sets a random seed based on the system and wall clock time
   implicit none

   integer , dimension(8)::val
   integer sze,i
   integer, allocatable, dimension(:):: pt
   logical ex

   integer :: clock_max
   integer :: clock_rate
   integer :: clock_reading

   call RANDOM_SEED(size=sze)
   allocate(pt(sze))

   inquire(file='seed.dat',exist=ex)
   if (ex) then
      open(10,file='seed.dat')
         read(10,*)pt
      close(10)
   else
      do i=1,sze
         call DATE_AND_TIME(values=val)
         call system_clock(clock_reading, clock_rate, clock_max)
         pt(i) = val(8)*val(7)+clock_reading
      enddo

      open(10,file='seed.dat')
         write(10,*)pt
      close(10)
   endif
   call RANDOM_SEED(put=pt)
   deallocate(pt)
end subroutine set_random_seed3
end module m_set_random_seed3
