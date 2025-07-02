module m_applyturbulence
contains
subroutine applyturbulence(f,turb_df,tau)
   use m_turbulenceforcing, only : iturb_pos, iturb_radius
   use m_readinfile, only : iforce
   use mod_dimensions
   use m_wtime
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)        ! distribution
   real, intent(in)    :: turb_df(nl,-ieps:ieps,ny,nz)      ! forcing distributions
   real, intent(in)    :: tau(nx,ny,nz)                     ! Tau
   integer ip,i,j,k,itr
   integer, parameter :: icpu=7
   call cpustart()
   ip=iturb_pos
   itr=iturb_radius
   if (iforce == 12 .or. iforce == 8) then
      do k=1,nz
      do j=1,ny
      do i=-iturb_radius,iturb_radius
         f(:,ip+i,j,k)=f(:,ip+i,j,k) + (1.0-0.5/tau(ip+i,j,k))*turb_df(:,i,j,k)
      enddo
      enddo
      enddo
   else
      f(1:nl,ip-itr:ip+itr,1:ny,1:nz) = f(1:nl,ip-itr:ip+itr,1:ny,1:nz) + turb_df(1:nl,-itr:itr,1:ny,1:nz)
   endif
   call cpufinish(icpu)
end subroutine
end module
