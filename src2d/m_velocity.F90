module m_velocity
! Computes density as a sum over particles with different velocities
contains
function velocity(f,rho,cs) result(vel)
   use mod_dimensions
   real,    intent(in) :: f(nx,ny,nl)
   real,    intent(in) :: rho(nx,ny)
   integer, intent(in) :: cs(nl)
   real vel(nx,ny)
   integer i,j,l

   vel=0.0
   do l=1,nl
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j) SHARED(l, vel, f, cs, rho)
   do j=1,ny
   do i=1,nx
      vel(i,j) =  vel(i,j) + f(i,j,l)*real(cs(l))/rho(i,j)
   enddo
   enddo
!$OMP END PARALLEL DO
   enddo

end function
end module
