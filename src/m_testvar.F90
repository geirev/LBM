module m_testvar
contains
subroutine testvar(variable,varname,ntot,it)
   implicit none
   integer, intent(in)          :: ntot
   real, intent(in)             :: variable(ntot)
   character(len=*), intent(in) :: varname
   real                         :: tmp(ntot)
   integer, intent(in)          :: it
#ifdef _CUDA
   attributes(device) :: variable
   attributes(device) :: tmp
#endif
   real(kind=4)                 :: variable_h(ntot)
   character(len=100) :: filename
   character(len=3)   :: cit
   logical ex
   integer iunit,i
   real fsum
   real fmax
   real eps
   real diff

   write(cit,'(i3.3)')it

   print '(a,a)','Testing (testvar):',trim(varname)

   filename='testvar_'//trim(varname)//cit//'.uf'

   eps = sqrt(tiny(1.0))


   inquire(file=trim(filename),exist=ex)

   if (ex) then
      open(newunit=iunit,file=trim(filename),form="unformatted", status='old')
         read(iunit)variable_h
      close(iunit)
      tmp=variable_h
#ifdef _CUDA
!$cuf kernel do(1) <<<*,*>>>
#endif
     fsum=0.0
     fmax=0.0
     do i=1,ntot
        diff=abs(variable(i)-tmp(i))
        if (diff > 1.0E-4) then
           print *,i,diff
        endif
        fsum=fsum+diff
        fmax=max(diff,fmax)
     enddo
     fsum=fsum/real(ntot)
     print *,'testvar: Mean abs error=',fsum,' max error=',fmax
   else
      print *,'testvar: saving file:',trim(filename)
      open(newunit=iunit,file=trim(filename),form="unformatted", status='replace')
         variable_h=variable
         write(iunit)variable_h
      close(iunit)
   endif

end subroutine
end module

