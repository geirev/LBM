module m_applyturbines
contains
subroutine applyturbines(f,df,tau)
   use mod_dimensions
   use m_readinfile, only : ipos,nturbines,iforce
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
   integer n,ip,i,j,k,ii
   integer, parameter :: icpu=8
   call cpustart()
   if (iforce == 12 .or. iforce == 8) then
      do n=1,nturbines
         ip=ipos(n)
#ifdef _CUDA
!$cuf kernel do(3) <<<*,*>>>
#endif
         do k=1,nz
         do j=1,ny
         do i=-ieps,ieps
            f(:,ip+i,j,k)=f(:,ip+i,j,k) + (1.0-0.5/tau(ip+i,j,k))*df(:,i,j,k,n)
         enddo
         enddo
         enddo
      enddo
   else
      do n=1,nturbines
         ip=ipos(n)
#ifdef _CUDA
!$cuf kernel do(2) <<<*,*>>>
#endif
         do k=1,nz
         do j=1,ny
         do ii=-ieps,ieps
            i=ip+ii
            !print '(4i4,10g13.4)',k,j,ii,i,df(1:10,ii,j,k,n)
            f(1:nl,i,j,k) = f(1:nl,i,j,k) + df(1:nl,ii,j,k,n)
         enddo
         enddo
         enddo
      enddo
   endif
   call cpufinish(icpu)
end subroutine
end module
