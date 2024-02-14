module m_density
! Computes density as a sum over particles with different velocities
contains
function density(f,blanking) result(dens)
   use mod_dimensions
   real,    intent(in) :: f(0:nx+1,ny,nz,nl)
   logical, intent(in) :: blanking(nx,ny,nz)
   real dens(nx,ny,nz)
   integer i,j,k,l
   dens(1:nx,:,:)=sum(f(1:nx,:,:,:),dim=4)
!   dens=0.0
!   do k=1,nz
!      do j=1,ny
!         do i=1,nx
!            if (.not.blanking(i,j,k)) dens(i,j,k)=dens(i,j,k)+f(i,j,k,l)
!         enddo
!      enddo
!    enddo
end function
end module
