module m_macrovars
! Computes density as a sum over particles with different velocities
contains
subroutine macrovars(rho,u,v,w,f,blanking)
   use mod_dimensions
   use m_wtime
   real,    intent(in)  :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(out) :: rho(nx,ny,nz)
   real,    intent(out) :: u(nx,ny,nz)
   real,    intent(out) :: v(nx,ny,nz)
   real,    intent(out) :: w(nx,ny,nz)
   logical, intent(in)  :: blanking(0:nx+1,0:ny+1,0:nz+1)
   integer i,j,k,l
   integer, parameter :: icpu=2
   if (nl /= 27) stop 'macrovars routine set up for D3Q27'

   call cpustart()


!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(u, v, w, rho, f, blanking)
   do k=1,nz
   do j=1,ny
   do i=1,nx
         rho(i,j,k)=f(1,i,j,k)
         do l = 2, nl
            rho(i,j,k)=rho(i,j,k)+f(l,i,j,k)
         enddo

         u(i,j,k) = f( 2,i,j,k)
         u(i,j,k) = u(i,j,k) - f( 3,i,j,k) &
                             + f( 8,i,j,k) &
                             - f( 9,i,j,k) &
                             + f(10,i,j,k) &
                             - f(11,i,j,k) &
                             - f(12,i,j,k) &
                             + f(13,i,j,k) &
                             - f(16,i,j,k) &
                             + f(17,i,j,k) &
                             - f(20,i,j,k) &
                             + f(21,i,j,k) &
                             - f(22,i,j,k) &
                             + f(23,i,j,k) &
                             + f(24,i,j,k) &
                             - f(25,i,j,k) &
                             - f(26,i,j,k) &
                             + f(27,i,j,k)
         u(i,j,k) =  u(i,j,k)/rho(i,j,k)

         v(i,j,k) = f( 4,i,j,k)
         v(i,j,k) = v(i,j,k) - f( 5,i,j,k) &
                             + f( 8,i,j,k) &
                             - f( 9,i,j,k) &
                             - f(10,i,j,k) &
                             + f(11,i,j,k) &
                             + f(14,i,j,k) &
                             - f(15,i,j,k) &
                             - f(18,i,j,k) &
                             + f(19,i,j,k) &
                             + f(20,i,j,k) &
                             - f(21,i,j,k) &
                             - f(22,i,j,k) &
                             + f(23,i,j,k) &
                             + f(24,i,j,k) &
                             - f(25,i,j,k) &
                             + f(26,i,j,k) &
                             - f(27,i,j,k)
         v(i,j,k) =  v(i,j,k)/rho(i,j,k)

         w(i,j,k) = -f( 6,i,j,k)
         w(i,j,k) = w(i,j,k) +f( 7,i,j,k) &
                             -f(12,i,j,k) &
                             +f(13,i,j,k) &
                             +f(14,i,j,k) &
                             -f(15,i,j,k) &
                             +f(16,i,j,k) &
                             -f(17,i,j,k) &
                             +f(18,i,j,k) &
                             -f(19,i,j,k) &
                             +f(20,i,j,k) &
                             -f(21,i,j,k) &
                             -f(22,i,j,k) &
                             +f(23,i,j,k) &
                             -f(24,i,j,k) &
                             +f(25,i,j,k) &
                             -f(26,i,j,k) &
                             +f(27,i,j,k)
         w(i,j,k) =  w(i,j,k)/rho(i,j,k)
!      endif
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

   call cpufinish(icpu)

end subroutine
end module
