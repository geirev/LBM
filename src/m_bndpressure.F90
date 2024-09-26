module m_bndpressure
! DOESN'T WORK
! Periodic boundary condition in i-direction with pressure difference
contains
subroutine bndpressure(f,rho,u,v,w)
   use mod_dimensions
   use m_fhrrscalar
   use m_readinfile
   implicit none
   real, intent(inout) :: f(0:nx+1,0:ny+1,0:nz+1,nl)
   real, intent(in)    :: rho(nx,ny,nz)
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)
   integer k,j,l

   real, allocatable :: fin(:,:,:)
   real, allocatable :: fou(:,:,:)

   real, allocatable :: fbin(:,:,:)
   real, allocatable :: fbou(:,:,:)

   allocate(fin(ny,nz,nl))
   allocate(fou(ny,nz,nl))
   allocate(fbin(ny,nz,nl))
   allocate(fbou(ny,nz,nl))

! Equilibrium distribution on boundaries
   do k=1,nz
   do j=1,ny
      fin(j,k,:)=fhrrscalar(rho(nx,j,k),u(nx,j,k),v(nx,j,k),w(nx,j,k))
      fou(j,k,:)=fhrrscalar(rho(1 ,j,k),u(1 ,j,k),v(1 ,j,k),w(1 ,j,k))
   enddo
   enddo

! New equilibrium distribution on boundaries
   do k=1,nz
   do j=1,ny
      fbin(j,k,:)=fhrrscalar(rho0+rhoa,u(nx,j,k),v(nx,j,k),w(nx,j,k))
      fbou(j,k,:)=fhrrscalar(rho0-rhoa,u(1 ,j,k),v(1 ,j,k),w(1 ,j,k))
   enddo
   enddo

! New boundary conditions
   do l=1,nl
   do k=1,nz
   do j=1,ny
      f(0   ,j,k,l)= fbin(j,k,l) + f(0   ,j,k,l) - fin(j,k,l)
      f(nx+1,j,k,l)= fbou(j,k,l) + f(nx+1,j,k,l) - fou(j,k,l)
   enddo
   enddo
   enddo
   deallocate(fin,fou,fbin,fbou)

end subroutine
end module
