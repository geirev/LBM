module m_turbines_apply
contains
subroutine turbines_apply(f,df,tau)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_turbines_init, only : ieps
   use m_readinfile   , only : ipos,nturbines,iforce
   use m_wtime
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)        ! distribution
   real, intent(in)    :: df(nl,-ieps:ieps,ny,nz,nturbines) ! forcing distributions
   real, intent(in)    :: tau(nx,ny,nz)                     ! Tau
#ifdef _CUDA
   attributes(device) :: f
   attributes(device) :: df
   attributes(device) :: tau
#endif
   integer n,ip,i,j,k,ii,l
   integer, parameter :: icpu=8
   real fac
   call cpustart()
   fac=1.0 ! for iforce=10
   do n=1,nturbines
      ip=ipos(n)
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
      do k=1,nz
      do j=1,ny
      do ii=-ieps,ieps
         i=ip+ii
         if (iforce == 8) fac=(1.0-0.5/tau(i,j,k))
         do l=1,nl
            f(l,i,j,k) = f(l,i,j,k) + fac*df(l,ii,j,k,n)
!            f(l,i,j,k) = f(l,i,j,k) + df(l,ii,j,k,n)
         enddo
      enddo
      enddo
      enddo
   enddo
   call cpufinish(icpu)
end subroutine
end module
