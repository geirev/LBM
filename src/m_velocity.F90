module m_velocity
! Computes density as a sum over particles with different velocities
contains
function velocity(f,rho,cs,blanking) result(vel)
   use mod_dimensions
   real,    intent(in) :: f(nx,ny,nz,nl)
   real,    intent(in) :: rho(nx,ny,nz)
   integer, intent(in) :: cs(nl)
   logical, intent(in) :: blanking(nx,ny,nz)
   real vel(nx,ny,nz)
   integer i,j,k,l

   vel=0.0
   do l=1,nl
!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(l, vel, f, cs, rho)
   do k=1,nz
   do j=1,ny
   do i=1,nx
      vel(i,j,k) =  vel(i,j,k) + f(i,j,k,l)*real(cs(l))/rho(i,j,k)
   enddo
   enddo
   enddo
!$OMP END PARALLEL DO
   enddo

!$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(i,j,k) SHARED(vel,blanking)
   do k=1,nz
      do j=1,ny
         do i=1,nx
            if (blanking(i,j,k)) then
               vel(i,j,k)=0.0
            endif
         enddo
      enddo
    enddo
!$OMP END PARALLEL DO

end function
end module
