module m_saverestart
contains
subroutine saverestart(it,f,theta,uu,vv,ww,rr,ibnd)
   use mod_dimensions
   integer, intent(in) :: it
   real,    intent(in) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in) :: theta
   real,    intent(in) :: uu(ny,nz,0:nrturb)
   real,    intent(in) :: vv(ny,nz,0:nrturb)
   real,    intent(in) :: ww(ny,nz,0:nrturb)
   real,    intent(in) :: rr(ny,nz,0:nrturb)
   integer, intent(in) :: ibnd

   character(len=6) cit
   integer :: irec

   write(cit,'(i6.6)')it
   print '(a,a)',' saverestart:',cit

   if (ibnd == 1) then
      inquire(iolength=irec)ny,nz,nrturb,uu,vv,ww,rr
      open(10,file='turbulence'//cit//'.uf',form="unformatted", access="direct", recl=irec)
         write(10,rec=1)ny,nz,nrturb,uu,vv,ww,rr
      close(10)
   endif

   open(10,file='theta'//cit//'.dat')
      write(10,*)theta
   close(10)

   inquire(iolength=irec) nx,ny,nz,nl,f
   open(10,file='restart'//cit//'.uf',form="unformatted", access="direct", recl=irec)
      write(10,rec=1)nx,ny,nz,nl,f
   close(10)

end subroutine
end module

