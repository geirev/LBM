module m_initialization
contains
function initialization(rho0,afac,ycos,zcos) result(f)
   use mod_dimensions
   use m_density
   real, intent(in)    :: rho0
   real, intent(in)    :: afac
   real, intent(in)    :: ycos
   real, intent(in)    :: zcos
   real f(nx,ny,nz,nl)
   real, allocatable :: rho(:,:,:)
   integer i,j,k
   real,    parameter :: pi   = 3.1415927  ! pi


   f=0.0

   ! intitialize to one + small perturbation
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

!   do k=1,nz
!      print '(a,i3,2f10.3)','z',k, afac*bcos*cos(2.0*pi*real(k-1)/real(nz)),f(1,1,k,2)
!   enddo
!   do j=1,ny
!      print '(a,i3,2f10.3)','y',j, afac*bcos*cos(2.0*pi*real(j-1)/real(ny)),f(1,j,1,2)
!   enddo

! Rescale to get constant reference density rho0
   allocate(rho(nx,ny,nz))
   rho=density(f)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      f(i,j,k,:)=f(i,j,k,:)*rho0/rho(i,j,k)
   enddo
   enddo
   enddo
   deallocate(rho)

end function
end module
