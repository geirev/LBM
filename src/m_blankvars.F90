module m_blankvars
contains
subroutine blankvars(u,v,w,vel,vortx,vorty,vortz,vort,blanking)
   use mod_dimensions
   implicit none
   real, intent(inout) :: u(nx,ny,nz)
   real, intent(inout) :: v(nx,ny,nz)
   real, intent(inout) :: w(nx,ny,nz)
   real, intent(inout) :: vel(nx,ny,nz)
   real, intent(inout) :: vortx(nx,ny,nz)
   real, intent(inout) :: vorty(nx,ny,nz)
   real, intent(inout) :: vortz(nx,ny,nz)
   real, intent(inout) :: vort(nx,ny,nz)
   logical, intent(in) :: blanking(nx,ny,nz)

   where (blanking)
      u=0.0
      v=0.0
      w=0.0
      vel=0.0
!      vortx=0.0
!      vorty=0.0
!      vortz=0.0
!      vort=0.0
   endwhere

end subroutine
end module
