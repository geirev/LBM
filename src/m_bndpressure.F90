module m_bndpressure
! DOESN'T WORK
! Periodic boundary condition in i-direction with pressure difference
contains
subroutine bndpressure(f,rho,u,v,w)
   use mod_dimensions
   use m_fequilscalar
   use m_readinfile
   implicit none
   real, intent(inout) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real, intent(in)    :: rho(nx,ny,nz)
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)
   integer k,j

   real, allocatable :: fin(:,:,:)
   real, allocatable :: fou(:,:,:)

   real, allocatable :: fbin(:,:,:)
   real, allocatable :: fbou(:,:,:)

   allocate( fin(nl,ny,nz))
   allocate( fou(nl,ny,nz))
   allocate(fbin(nl,ny,nz))
   allocate(fbou(nl,ny,nz))

! Equilibrium distribution on boundaries
   do k=1,nz
   do j=1,ny
      fin(:,j,k)=fequilscalar(rho(nx,j,k),u(nx,j,k),v(nx,j,k),w(nx,j,k))
      fou(:,j,k)=fequilscalar(rho(1 ,j,k),u(1 ,j,k),v(1 ,j,k),w(1 ,j,k))
   enddo
   enddo

! New equilibrium distribution on boundaries
   do k=1,nz
   do j=1,ny
      fbin(:,j,k)=fequilscalar(rho0+rhoa,u(nx,j,k),v(nx,j,k),w(nx,j,k))
      fbou(:,j,k)=fequilscalar(rho0-rhoa,u(1 ,j,k),v(1 ,j,k),w(1 ,j,k))
   enddo
   enddo

! New boundary conditions
   do k=1,nz
   do j=1,ny
      f(:,0   ,j,k)= fbin(:,j,k) + f(:,0   ,j,k) - fin(:,j,k)
      f(:,nx+1,j,k)= fbou(:,j,k) + f(:,nx+1,j,k) - fou(:,j,k)
   enddo
   enddo
   deallocate(fin,fou,fbin,fbou)

end subroutine
end module
