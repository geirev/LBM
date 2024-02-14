module m_bndpressure
! DOESN'T WORK
! Periodic boundary condition in i-direction with pressure difference
contains
subroutine bndpressure(f,rho,u,v,w)
   use mod_dimensions
   use m_feqscalar
   implicit none
   real, intent(inout) :: f(0:nx+1,ny,nz,nl)
   real, intent(in)    :: rho(nx,ny,nz)
   real, intent(in)    :: u(nx,ny,nz)
   real, intent(in)    :: v(nx,ny,nz)
   real, intent(in)    :: w(nx,ny,nz)
   integer k,j,l

   real fin(ny,nz,nl)
   real fou(ny,nz,nl)

   real fbin(ny,nz,nl)
   real fbou(ny,nz,nl)

! Equilibrium distribution on boundaries
   do k=1,nz
   do j=1,ny
      call feqscalar(fin(j,k,:),rho(nx,j,k),u(nx,j,k),v(nx,j,k),w(nx,j,k))
      call feqscalar(fou(j,k,:),rho(1 ,j,k),u(1 ,j,k),v(1 ,j,k),w(1 ,j,k))
   enddo
   enddo

! New equilibrium distribution on boundaries
   do k=1,nz
   do j=1,ny
      call feqscalar(fbin(j,k,:),101.0,u(nx,j,k),v(nx,j,k),w(nx,j,k))
      call feqscalar(fbou(j,k,:), 99.0,u(1 ,j,k),v(1 ,j,k),w(1 ,j,k))
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

end subroutine
end module
