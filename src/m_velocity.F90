module m_velocity
! Computes density as a sum over particles with different velocities
contains
function velocity(f,rho,cs,blanking) result(vel)
   use mod_dimensions
   use m_wtime
   real,    intent(in) :: f(nl,0:nx+1,0:ny+1,0:nz+1)
   real,    intent(in) :: rho(nx,ny,nz)
   integer, intent(in) :: cs(nl)
   logical, intent(in) :: blanking(0:nx+1,0:ny+1,0:nz+1)
   real vel(nx,ny,nz)
   integer i,j,k,l
   integer, parameter :: icpu=2
   call cpustart()


   vel=0.0
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k,l) SHARED(vel, f, cs)
   do k=1,nz
   do j=1,ny
   do i=1,nx
   do l=1,nl
      vel(i,j,k) =  vel(i,j,k) + f(l,i,j,k)*real(cs(l))
   enddo
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(vel, cs, rho, blanking)
   do k=1,nz
   do j=1,ny
   do i=1,nx
     !if (blanking(i,j,k)) then
     !   vel(i,j,k)=0.0
     !else
         vel(i,j,k) =  vel(i,j,k)/rho(i,j,k)
     !endif
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
   call cpufinish(icpu)

end function
end module
