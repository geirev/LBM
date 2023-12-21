module   m_initialization
contains
function initialization(rho0,afac,bcos) result(f)
   use mod_dimensions
   use m_density
   real, intent(in)    :: rho0
   real, intent(in)    :: afac
   real, intent(in)    :: bcos
   real f(nx,ny,nl)
   real, allocatable :: rho(:,:)
   integer i,j
   real,    parameter :: pi   = 3.1415927  ! pi

!   print '(a,3i4)','Initialization',nx,ny,nl

   f=0.0

   ! intitialize to one + small perturbation
   call random_number(f)
   f=1.0 + 0.01*f

   ! Adding horizontal velocity component in x direction with some y-variation
   do j=1,ny
      f(:,j,4) =  f(:,j,4) + afac*(1.0 + bcos*cos(2.0*pi*real(j-1)/ny))
   enddo


! Rescale to get constant reference density rho0
   allocate(rho(nx,ny))
   rho=density(f)
   do j=1,ny
   do i=1,nx
      f(i,j,:)=f(i,j,:)*rho0/rho(i,j)
   enddo
   enddo
   deallocate(rho)

end function
end module
