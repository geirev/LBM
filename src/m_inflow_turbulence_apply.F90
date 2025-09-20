module m_inflow_turbulence_apply
! only set up for Kupershtokh forcing
contains
subroutine inflow_turbulence_apply(f,turbulence_df)
   use m_inflow_turbulence_init, only : iturb_pos, iturb_radius
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_wtime
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)   ! distribution
   real, intent(in)    :: turbulence_df(nl,ny,nz)      ! forcing distributions
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: turbulence_df
#endif
   integer ip,j,k,l
   integer, parameter :: icpu=9
   call cpustart()
   ip=iturb_pos
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
   do k=1,nz
   do j=1,ny
   do l=1,nl
      f(l,ip,j,k) = f(l,ip,j,k) + turbulence_df(l,j,k)
   enddo
   enddo
   enddo
   call cpufinish(icpu)
end subroutine
end module
