module m_density
! Computes density as a sum over particles with different velocities
contains
function density(f,blanking) result(dens)
   use mod_dimensions
   use m_wtime
   real,    intent(in) :: f(0:nx+1,0:ny+1,0:nz+1,nl)
   logical, intent(in) :: blanking(nx,ny,nz)
   real dens(nx,ny,nz)
   integer i,j,k,l
   integer, parameter :: icpu=3
   call cpustart()

!   dens(1:nx,1:ny,1:nz)=sum(f(1:nx,1:ny,1:nz,:),dim=4)

   dens(:,:,:)=0.0
   do l = 1, nl
!$OMP PARALLEL DO PRIVATE(i,j,k) SHARED(f, dens, l)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      dens(i,j,k)=dens(i,j,k)+f(i,j,k,l)
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
   enddo

   call cpufinish(3)
end function
end module
