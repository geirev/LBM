module m_seedmanagement
contains
subroutine seedmanagement(nt0)
   use m_set_random_seed3
   integer, intent(in) :: nt0
   logical ex

   inquire(file='seed.orig',exist=ex)
   if (nt0 == 0) then
      if (ex) then
         print '(a)','Copying seed.orig to seed.dat and loading seed from seed.dat'
         call system('cp seed.orig seed.dat')
         call set_random_seed3
      else
         print '(a)','Generating new seed.dat and copying it to seed.orig'
         call set_random_seed3
         call system('cp seed.dat seed.orig')
      endif
   endif

   if (nt0 > 0) then
      inquire(file='seed.dat',exist=ex)
      if (ex) then
         print '(a)','Removing old seed.dat and setting new seed stored in seed.dat'
         call system('rm seed.dat')
         call set_random_seed3
      endif
   endif

end subroutine
end module
