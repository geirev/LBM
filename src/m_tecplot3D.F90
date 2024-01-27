module m_tecplot3D
contains
subroutine tecplot3D(fname,rho,u,v,w,vel,vortx,vorty,vortz,vort,blanking)
   use mod_dimensions
   implicit none

   character(len=*), intent(in)     :: fname
   logical, intent(in)  :: blanking(nx,ny,nz)
   real, intent(in)     :: rho(nx,ny,nz)
   real, intent(in)     :: u(nx,ny,nz)
   real, intent(in)     :: v(nx,ny,nz)
   real, intent(in)     :: w(nx,ny,nz)
   real, intent(in)     :: vel(nx,ny,nz)
   real, intent(in)     :: vortx(nx,ny,nz)
   real, intent(in)     :: vorty(nx,ny,nz)
   real, intent(in)     :: vortz(nx,ny,nz)
   real, intent(in)     :: vort(nx,ny,nz)
   integer iblanking(nx,ny,nz)
   integer i,j,k

   iblanking=0
   where (blanking) iblanking=1

   open(10,file=trim(fname))
   write(10,'(3a)') 'TITLE = "LBM3D output ',trim(fname),'"'
   write(10,'(a)',advance='no')'VARIABLES = "I-index" "J-index" "K-index"'
   write(10,'(a)')' "blanking" "rho" "u" "v" "w" "vel" "vortx" "vorty" "vortz" "vort"'

   write(10,'(a,a,3(a,i5))')'ZONE T="',trim(fname),'"  F=BLOCK, I=',nx,' ,J=',ny,' ,K=',nz

   write(10,'(30i5)')(((i,i=1,nx),j=1,ny),k=1,nz)
   write(10,'(30i5)')(((j,i=1,nx),j=1,ny),k=1,nz)
   write(10,'(30i5)')(((k,i=1,nx),j=1,ny),k=1,nz)
   write(10,'(30i2)')iblanking
   write(10,'(10(1x,g15.7))')rho
   write(10,'(10(1x,g15.7))')u
   write(10,'(10(1x,g15.7))')v
   write(10,'(10(1x,g15.7))')w
   write(10,'(10(1x,g15.7))')vel
   write(10,'(10(1x,g15.7))')vortx
   write(10,'(10(1x,g15.7))')vorty
   write(10,'(10(1x,g15.7))')vortz
   write(10,'(10(1x,g15.7))')vort
   close(10)

end subroutine
end module

