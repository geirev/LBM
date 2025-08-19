module m_applyturbulence
contains
subroutine applyturbulence(f,turbulence_df,tau)
   use m_turbulenceforcing, only : iturb_pos, iturb_radius
   use m_readinfile, only : iforce
   use mod_dimensions
   use m_wtime
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)        ! distribution
   real, intent(in)    :: turbulence_df(nl,ny,nz)      ! forcing distributions
   real, intent(in)    :: tau(nx,ny,nz)                     ! Tau
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: turbulence_df
   attributes(device) :: tau
#endif
   real tmp
   integer ip,i,j,k,l
   integer, parameter :: icpu=9
   call cpustart()
   ip=iturb_pos
   if (iforce == 12 .or. iforce == 8) then
      stop 'only works for iforce=10'
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,nz
      do j=1,ny
      do i=-iturb_radius,iturb_radius
         f(:,ip+i,j,k)=f(:,ip+i,j,k) + (1.0-0.5/tau(ip+i,j,k))*turbulence_df(:,j,k)
      enddo
      enddo
      enddo
   else
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=1,nz
      do j=1,ny
      do l=1,nl
         tmp=turbulence_df(l,j,k)
         f(l,ip,j,k) = f(l,ip,j,k) + tmp
      enddo
      enddo
      enddo
   endif
   call cpufinish(icpu)
end subroutine
end module
