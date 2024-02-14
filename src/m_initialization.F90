module m_initialization
contains
subroutine initialization(f,blanking)
   use mod_dimensions
   use m_density
   use m_readinfile
   real,    intent(inout)   :: f(0:nx+1,ny,nz,nl)
   logical, intent(inout) :: blanking(nx,ny,nz)
   real,    allocatable   :: rho(:,:,:)
   real,    parameter     :: pi   = 3.1415927  ! pi
   integer i,j,k



   ! intitialize to feq + small perturbation
   call random_number(f)
   f=1.0 + 0.01*f

   ! Adding horizontal velocity component in x direction with some y-variation
   do k=1,nz
   do j=1,ny
      f(:,j,k,2) =  f(:,j,k,2) +                                    &
                  afac*(1.0 + ycos*cos(2.0*pi*real(j-1)/real(ny))   &
                            + zcos*cos(2.0*pi*real(k-1)/real(nz)) )
   enddo
   enddo

! Rescale to get constant reference density rho0
   allocate(rho(nx,ny,nz))
   rho=density(f,blanking)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      f(i,j,k,:)=f(i,j,k,:)*rho0/rho(i,j,k)
   enddo
   enddo
   enddo
   deallocate(rho)

end subroutine
end module
