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
   integer n,ip,i,j,k
   integer, parameter :: icpu=7
   call cpustart()
   if (iforce == 12 .or. iforce == 8) then
      do n=1,nturbines
         ip=ipos(n)
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
         f(1:nl,ip-ieps:ip+ieps,1:ny,1:nz) = f(1:nl,ip-ieps:ip+ieps,1:ny,1:nz) + df(1:nl,-ieps:ieps,1:ny,1:nz,n)
      enddo
   endif
   call cpufinish(icpu)
end subroutine
end module
