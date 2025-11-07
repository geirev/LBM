module m_seedmanagement
contains
subroutine seedmanagement(nt0,irank)
   use m_set_random_seed3
   implicit none
   integer, intent(in) :: nt0
   integer, intent(in) :: irank
   logical ex
   character(len=4) :: crank
   character(len=100) :: fnamed
   character(len=100) :: fnameo

   write(crank,'(i4.4)')irank
   fnameo='seed_'//crank//'.orig'
   fnamed='seed_'//crank//'.dat'
   inquire(file=trim(fnameo),exist=ex)
   if (nt0 == 0) then
      if (ex) then
         print '(a)','Copying seed.orig to seed.dat and loading seed from seed.dat'
         call system('cp '//trim(fnameo)//' '//trim(fnamed))
         call set_random_seed3(irank)
      else
         print '(a)','Generating new seed.dat and copying it to seed.orig'
         call set_random_seed3(irank)
         call system('cp '//trim(fnamed)//' '//trim(fnameo))
      endif
   endif

   if (nt0 > 0) then
      inquire(file='seed.dat',exist=ex)
      if (ex) then
         print '(a)','Removing old seed.dat and setting new seed stored in seed.dat'
         call system('rm '//trim(fnamed))
         call set_random_seed3(irank)
      endif
   endif

end subroutine
end module
