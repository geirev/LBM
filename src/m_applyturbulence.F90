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
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: turb_df
   attributes(device) :: tau
#endif
   integer ip,i,j,k,l,itr,ii
   integer, parameter :: icpu=7
   call cpustart()
   ip=iturb_pos
   itr=iturb_radius
   if (iforce == 12 .or. iforce == 8) then
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,nz
      do j=1,ny
      do i=-iturb_radius,iturb_radius
         f(:,ip+i,j,k)=f(:,ip+i,j,k) + (1.0-0.5/tau(ip+i,j,k))*turb_df(:,i,j,k)
      enddo
      enddo
      enddo
   else
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
      do k=1,nz
      do j=1,ny
      do ii=-itr,itr
      i=ip+ii
      do l=1,nl
         f(l,i,j,k) = f(l,i,j,k) + turb_df(l,ii,j,k)
      enddo
      enddo
      enddo
      enddo
   endif
   call cpufinish(icpu)
end subroutine
end module
