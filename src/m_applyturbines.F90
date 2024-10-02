module m_applyturbines
contains
subroutine applyturbines(f,df,tau)
   use mod_dimensions
   use m_readinfile, only : ipos,nturbines,iforce
   use m_wtime
   implicit none
   real, intent(inout) :: f(0:nx+1,0:ny+1,0:nz+1,nl)        ! distribution
   real, intent(in)    :: df(-ieps:ieps,ny,nz,nl,nturbines) ! forcing distributions
   real, intent(in)    :: tau(nx,ny,nz)                     ! Tau
   integer n,ip,i,j,k
   integer, parameter :: icpu=7
   call cpustart()
   if (iforce == 12) then
      do n=1,nturbines
         ip=ipos(n)
         do k=1,nz
         do j=1,ny
         do i=-ieps,ieps
            f(ip+i,j,k,:)=f(ip+i,j,k,:) + (1.0-0.5/tau(ip+i,j,k))*df(i,j,k,:,n)
         enddo
         enddo
         enddo
      enddo
   else
      do n=1,nturbines
         ip=ipos(n)
         f(ip-ieps:ip+ieps,1:ny,1:nz,1:nl) = f(ip-ieps:ip+ieps,1:ny,1:nz,1:nl) + df(-ieps:ieps,1:ny,1:nz,1:nl,n)
      enddo
   endif
   call cpufinish(icpu)
end subroutine
end module
