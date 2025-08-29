module m_density
! Computes density as a sum over particles with different velocities
contains
function density(f,blanking) result(dens)
   use mod_dimensions
   use mod_D3Q27setup, only : nl
   use m_wtime
   implicit none
   real,    intent(in) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   logical, intent(in) :: blanking(0:nx+1,0:ny+1,0:nz+1)
   real dens(nx,ny,nz)
   integer i,j,k,l
   integer, parameter :: icpu=1
   call cpustart()


   dens(:,:,:)=0.0
!$OMP PARALLEL DO PRIVATE(i,j,k,l) SHARED(f, dens)
   do k=1,nz
   do j=1,ny
   do i=1,nx
   do l = 1, nl
      dens(i,j,k)=dens(i,j,k)+f(l,i,j,k)
   enddo
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

   call cpufinish(icpu)
end function
end module
